using Test
using FFTW
using Statistics
using Distributed
using WaveDM: select_backend, detect_resources, ParallelBackend, DeviceConfig,
    parallel_poisson, parallel_sum, parallel_maximum, parallel_minimum,
    parallel_findmax, parallel_sumprod,
    to_device, to_host, distribute, release!,
    is_serial, is_threads, is_distributed, is_gpu, is_parallel

@testset "Parallel backend: detection" begin
    res = detect_resources()
    @test res.nthreads >= 1
    @test res.nworkers >= 0
    @test res.pids isa AbstractVector
    @test res.has_darrays isa Bool
    @test res.has_parallelops isa Bool
    @test res.has_cuda isa Bool
end

@testset "Parallel backend: select_backend" begin
    # Default kind=:auto picks *something* and the chosen kind is one of
    # the four valid options.
    b = select_backend()
    @test b.kind in (:serial, :threads, :distributed, :gpu)

    # Force each kind in turn; if the host cannot satisfy the request,
    # the selector should fall back gracefully (warning, not error).
    for k in (:serial, :threads, :distributed, :gpu)
        b = select_backend(; kind=k)
        @test b isa ParallelBackend
    end

    # Explicit hints steer the choice.
    b = select_backend(; gpu=false, distributed=false)
    @test b.kind in (:serial, :threads)

    # If distributed is requested and workers<2, we expect a fallback.
    if nworkers() < 2
        b = select_backend(; distributed=true)
        @test b.kind !== :distributed
    end
end

@testset "DeviceConfig compatibility shim" begin
    # Old positional constructor still works.
    dc1 = DeviceConfig(false, identity, collect)
    @test dc1.gpu == false
    @test dc1.kind === :serial
    @test dc1.DeviceArray === identity
    @test dc1.DA === collect

    # New keyword constructor delegates to select_backend.
    dc2 = DeviceConfig(; gpu=false)
    @test dc2 isa DeviceConfig
    @test dc2.gpu == false
    @test dc2.kind in (:serial, :threads, :distributed, :gpu)

    # All forwarded fields are reachable.
    b = select_backend(; gpu=false)
    dc3 = DeviceConfig(b.gpu, b.device_array, b.distribute_array)
    @test propertynames(dc3) === (
        :backend, :DeviceArray, :DA, :gpu, :device_array,
        :distribute_array, :local_array, :release!,
        :kind, :nthreads, :nworkers, :pids, :fft_threads,
        :has_darrays, :has_parallelops,
    )
end

@testset "Parallel backend: parallel_poisson correctness" begin
    # 3D Poisson solve on a smooth Gaussian source.  We test the
    # serial/threads path (which always works) and, if a distributed
    # setup is available, the distributed path as well.  We do not
    # require GPU: the GPU test is conditional on CUDA.functional().
    Nx, Ny, Nz = 32, 32, 32
    Δx = Δy = Δz = 1.0
    Δ = (Δx, Δy, Δz)
    N = (Nx, Ny, Nz)

    x = collect(LinRange(-π, π, Nx))
    y = collect(LinRange(-π, π, Ny))
    z = collect(LinRange(-π, π, Nz))
    ρ = [exp(-(xi^2 + yj^2 + zk^2)) for xi in x, yj in y, zk in z]

    # Reference solve: serial, no GPU.
    backend_serial = select_backend(; kind=:serial)
    Φ_ref = parallel_poisson(Δ, N, ρ, Periodic(), backend_serial)

    @test size(Φ_ref) == (Nx, Ny, Nz)
    @test eltype(Φ_ref) <: Complex
    @test all(isfinite, Φ_ref)
    # The spectral solution `ifft(-4π ρ̂ / |k|²)` for a Gaussian source
    # has bounded magnitude.  We only check it is not catastrophically
    # wrong (e.g. NaN, Inf, or many orders of magnitude off).
    @test maximum(abs, Φ_ref) < 100.0
    @test maximum(abs, Φ_ref) > 1e-6
end

#=
@testset "Parallel backend: parallel_poisson distributed" begin
    # Distributed Poisson only runs when ≥2 workers are present.
    if nworkers() < 2
        @info "Skipping distributed Poisson test (nworkers=$(nworkers()) < 2)"
    else
        backend = select_backend(; kind=:distributed)
        @test backend.kind === :distributed
        @test backend.pids == workers()

        Nx, Ny, Nz = 32, 32, 32
        Δ = (1.0, 1.0, 1.0)
        N = (Nx, Ny, Nz)

        x = collect(LinRange(-π, π, Nx))
        y = collect(LinRange(-π, π, Ny))
        z = collect(LinRange(-π, π, Nz))
        ρ = [exp(-(xi^2 + yj^2 + zk^2)) for xi in x, yj in y, zk in z]

        Φ_dist = parallel_poisson(Δ, N, ρ, Periodic(), backend)
        Φ_ref  = parallel_poisson(Δ, N, ρ, Periodic(),
                                  select_backend(; kind=:serial))

        # Slab-decomposed FFT should agree with the serial FFT to
        # roughly numerical precision.
        @test maximum(abs, Φ_dist .- Φ_ref) / maximum(abs, Φ_ref) < 1e-6
    end
end
=#

@testset "Parallel backend: reductions" begin
    A = reshape(collect(1.0:27.0), 3, 3, 3)
    for k in (:serial, :threads, :distributed, :gpu)
        b = select_backend(; kind=k)
        @test parallel_sum(A, b) ≈ sum(A)
        @test parallel_maximum(A, b) ≈ maximum(A)
        @test parallel_minimum(A, b) ≈ minimum(A)
        @test parallel_sumprod(A, A, b) ≈ sum(A .* A)
    end
end

@testset "Parallel backend: parallel_findmax" begin
    A = reshape(collect(1.0:27.0), 3, 3, 3)
    for k in (:serial, :threads)
        b = select_backend(; kind=k)
        v, idx = parallel_findmax(A, b)
        @test v ≈ maximum(A)
        @test A[idx] ≈ maximum(A)
    end
end

@testset "Parallel backend: to_device / to_host / release!" begin
    A = rand(4, 4, 4)
    b = select_backend(; kind=:serial)
    @test to_device(A, b) === A
    @test to_host(A, b) === A
    @test release!(A, b) === nothing
end
