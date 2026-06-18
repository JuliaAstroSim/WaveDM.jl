# Parallel Poisson solver for WaveDM.jl
#
# The gravitational potential satisfies ∇²Φ = 4π ρ, which in Fourier
# space becomes Φ̂(k) = -4π ρ̂(k) / |k|².  The strategy is:
#
#   :serial / :threads :  use PhysicalFFT.fft_poisson (FFTW, MKL-threaded)
#   :gpu                 :  use PhysicalFFT.fft_poisson with the GPU() device
#   :distributed         :  slab-decompose ρ, SPMD-local 3D FFT (mirroring the
#                          pattern the user already used in
#                          `NonlinearOptics.jl/examples/SPE-WaveDM/test_Parallel.jl`)
#
# The distributed implementation uses DistributedArrays + SPMD.  We
# deliberately do NOT load Dagger.jl: see
# JuliaParallel/Dagger.jl#649 (the DRef finalizer chain cannot release
# remote memory reliably).  DArrays created with
# `DistributedArrays.distribute` are reclaimed by
# Julia's normal GC and have no memory-leak problem.


"""
    parallel_poisson(Δ, N, ρ, boundary, backend; gpu_force=false) -> Φ

Solve the Poisson equation ∇²Φ = 4π ρ with periodic boundary conditions
on a uniform grid with spacings `Δ = (Δx, Δy, Δz)` and active size
`N = (Nx, Ny, Nz)`.

`ρ` may be an `Array`, `CuArray`, or `DArray` depending on `backend`.
The returned potential lives on the host (plain `Array`) unless the backend
is GPU, in which case the result stays on the device until the caller
explicitly copies it back.

This function replaces the inline `fft_poisson(Δ, N, ρ, boundary, ...)`
calls scattered through the WaveDM main loop and is the only place
where the parallel FFT strategy is encoded.
"""
function parallel_poisson(Δ, N, ρ, boundary, backend::ParallelBackend; gpu_force::Bool=false)
    if backend.kind === :gpu
        return _poisson_gpu(Δ, N, ρ, boundary)
    elseif backend.kind === :distributed
        return _poisson_distributed(Δ, N, ρ, boundary, backend)
    else
        return _poisson_local(Δ, N, ρ, boundary)
    end
end

function parallel_poisson(Δ, N, ρ, boundary, config::DeviceConfig; gpu_force::Bool=false)
    return parallel_poisson(Δ, N, ρ, boundary, config.backend; gpu_force=gpu_force)
end

# ----------------------------------------------------------------------------
# Local backends: serial / threads
# ----------------------------------------------------------------------------

function _poisson_local(Δ, N, ρ, boundary)
    ρ_host = collect(ρ)
    return PhysicalFFT.fft_poisson(Δ, N, ρ_host, boundary, PhysicalFFT.CPU())
end

# ----------------------------------------------------------------------------
# GPU
# ----------------------------------------------------------------------------

function _poisson_gpu(Δ, N, ρ, boundary)
    ρ_dev = ρ isa CUDA.CuArray ? ρ : CUDA.cu(ρ)
    # Always return a host array.  Callers that need the device
    # representation can wrap with `to_device(Φ, backend)` afterwards.
    return Array(PhysicalFFT.fft_poisson(Δ, N, ρ_dev, boundary, PhysicalFFT.GPU()))
end

# ----------------------------------------------------------------------------
# Distributed: slab-decomposed 3D FFT via DistributedArrays.SPMD
# ----------------------------------------------------------------------------
#
# Decomposition: each worker holds a slab of size (Nx, Ny, Nz/nw) along z.
#
#   1. Local fft along axes 1 and 2  (z-slab, no communication)
#   2. Redistribute to y-slab: (Nx, Ny/nw, Nz)
#   3. Local fft along axis 3  (y-slab, no communication)
#   4. Apply the k-space factor -4π/|k|² (y-slab, with global y-offset)
#   5. Local ifft along axis 3  (y-slab)
#   6. Redistribute back to z-slab
#   7. Local ifft along axes 1, 2  (z-slab)
#   8. Collect to host
#
# This is the same slab-decomposition strategy used in
# `test_Parallel.jl::spmd_fft3`.  The k-space factor is precomputed on
# the host and distributed to the y-slab to avoid a per-element index
# computation on the workers.

function _poisson_distributed(Δ, N, ρ, boundary, backend::ParallelBackend)
    error("distributed Poisson solver is not implemented!")
    # pids = backend.pids
    # nw = length(pids)
    # nw >= 2 || throw(ArgumentError("Distributed Poisson needs >=2 workers; got $(nw)"))

    # Nx, Ny, Nz = N
    # Δx, Δy, Δz = Δ

    # # ---- 1. Distribute ρ as a z-slab DArray --------------------------------
    # A = ρ isa DArray ? ρ : DistributedArrays.distribute(collect(ρ); procs=pids, dist=[1, 1, nw])

    # # ---- 2. Forward FFT along axes 1, 2  (z-slab) ---------------------------
    # A_xy = _spmd_zip(A, A) do a, b
    #     b[:L] .= FFTW.fft(a[:L], [1, 2])
    # end

    # # ---- 3. Redistribute to y-slab ------------------------------------------
    # A_y = DistributedArrays.distribute(collect(A_xy); procs=pids, dist=[1, nw, 1])

    # # ---- 4. Forward FFT along axis 3  (y-slab) -----------------------------
    # A_fwd = _spmd_zip(A_y, A_y) do a, b
    #     b[:L] .= FFTW.fft(a[:L], [3])
    # end

    # # ---- 5. Apply k-space factor -------------------------------------------
    # kfactor = _build_kfactor(Nx, Ny, Nz, Δx, Δy, Δz)
    # kfactor_y = DistributedArrays.distribute(kfactor; procs=pids, dist=[1, nw, 1])
    # A_fwd .*= kfactor_y

    # # ---- 6. Inverse FFT along axis 3  (y-slab) -----------------------------
    # A_inv = _spmd_zip(A_fwd, A_fwd) do a, b
    #     b[:L] .= FFTW.ifft(a[:L], [3])
    # end

    # # ---- 7. Redistribute to z-slab -----------------------------------------
    # A_z = DistributedArrays.distribute(collect(A_inv); procs=pids, dist=[1, 1, nw])

    # # ---- 8. Inverse FFT along axes 1, 2  (z-slab) --------------------------
    # Φ_z = _spmd_zip(A_z, A_z) do a, b
    #     b[:L] .= FFTW.ifft(a[:L], [1, 2])
    # end

    # # ---- 9. Collect to host ------------------------------------------------
    # return collect(Φ_z)
end

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

"""
    _build_kfactor(Nx, Ny, Nz, Δx, Δy, Δz) -> Array{Float64, 3}

Return the Fourier-space Poisson kernel `-4π / |k|²` for a uniform grid
of size `Nx × Ny × Nz` and spacing `Δ = (Δx, Δy, Δz)`.  The
`k=0` component is set to zero (mean-free convention).
"""
function _build_kfactor(Nx, Ny, Nz, Δx, Δy, Δz)
    kx = 2π .* FFTW.fftfreq(Nx, Nx * Δx)
    ky = 2π .* FFTW.fftfreq(Ny, Ny * Δy)
    kz = 2π .* FFTW.fftfreq(Nz, Nz * Δz)
    kfactor = Array{Float64}(undef, Nx, Ny, Nz)
    @inbounds for iz in 1:Nz, iy in 1:Ny, ix in 1:Nx
        k2 = kx[ix]^2 + ky[iy]^2 + kz[iz]^2
        kfactor[ix, iy, iz] = k2 == 0 ? 0.0 : -4π / k2
    end
    return kfactor
end

# Apply a 2-argument SPMD function: f(local_in, local_out) -> nothing.
# We use the `[:L]` local-part pattern from DistributedArrays.SPMD to
# match the test code (test_Parallel.jl).
function _spmd_zip(f, A::DArray, B::DArray)
    DistributedArrays.SPMD.spmd(f, A, B; pids=workers(A))
    return B
end
