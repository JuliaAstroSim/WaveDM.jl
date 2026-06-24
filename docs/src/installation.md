# Installation

## Prerequisites

WaveDM.jl requires **Julia 1.11 or later**. Download Julia from the [official website](https://julialang.org/downloads/) and verify the installation:

```bash
julia --version
```

### Optional: GPU support

If you plan to use GPU acceleration you will need:

- An **NVIDIA GPU** with CUDA support.
- The **CUDA Toolkit** (version 12.x recommended) and a working `nvidia-smi` on `PATH`.
- On Windows, use the latest available NVIDIA driver that matches the CUDA toolkit bundled with `CUDA.jl`.

## Installation

### From the Julia registry (recommended)

Open a Julia REPL, press `]` to enter package mode, then run:

```julia
pkg> add WaveDM
```

### From GitHub (development version)

```julia
pkg> add https://github.com/JuliaAstroSim/WaveDM.jl
```

### From a local clone (developer workflow)

```bash
git clone https://github.com/JuliaAstroSim/WaveDM.jl.git
cd WaveDM.jl
julia --project -e 'using Pkg; Pkg.develop(PackageSpec(path="."))'
```

## Required dependencies

The following packages are required at runtime. They are installed automatically by `Pkg.add WaveDM`, so no extra step is needed for a CPU-only installation.

| Package | Role |
| --- | --- |
| `AstroIC` | Initial condition generation |
| `AstroNbodySim` | N-body gravity and the parallel runtime glue |
| `AstroPlot` | Analysis and post-processing visualization |
| `AstroSimBase` | Common type and unit definitions |
| `PhysicalParticles` | N-body particle types and vector algebra |
| `PhysicalFFT` | Spectral primitives for the SPE solver |
| `PhysicalFDM` | Finite-difference operators (gradient, boundary) |
| `PhysicalMeshes` | Mesh bookkeeping used by the FFT Poisson solver |
| `GalacticDynamics` | Density profiles and orbit calculations |
| `FFTW` | CPU FFT backend |
| `Distributed` | Distributed-memory parallelism |
| `DataFrames`, `CSV`, `JLD2` | Data I/O |
| `Unitful`, `UnitfulAstro` | Unit handling |
| `LsqFit`, `Optim`, `Dierckx` | Profile fitting |
| `Observables`, `StatsBase`, `SpecialFunctions` | Numerical utilities |
| `GLMakie`, `CairoMakie`, `UnicodePlots` | Visualization backends |

## Optional dependencies

| Package | Required for |
| --- | --- |
| `CUDA` | GPU acceleration (`gpu = true`) |
| `Sundials` | High-order N-body integrators via `AstroNbodySim` |
| `Dagger` | Distributed task graph (advanced) |
| `Tullio`, `Zygote` | Macro-based and AD-accelerated kernels (advanced) |
| `PaddedViews`, `Strided` | Strided-array optimizations for large grids |
| `OffsetArrays` | Offset indexing for ghost-zone arrays |

## Setting up GPU support

### 1. Install the CUDA toolkit

Follow the [`CUDA.jl` installation guide](https://cuda.juliagpu.org/stable/installation/overview/) for your platform. On Windows:

1. Install the latest NVIDIA Studio / Game Ready driver.
2. Verify that `nvidia-smi` runs and shows your GPU.
3. In the Julia REPL, run `using CUDA; CUDA.functional()?true` — it should print `true`.

### 2. Verify in Julia

```julia
using CUDA

if CUDA.has_cuda()
    println("CUDA is available")
    println("GPU:           ", CUDA.name(CUDA.device()))
    println("Capability:    ", CUDA.capability(CUDA.device()))
    println("Total memory:  ", CUDA.totalmemory(CUDA.device()) ÷ 1024^3, " GiB")
else
    println("CUDA is not available")
end
```

### 3. Smoke test WaveDM.jl on the GPU

```julia
using WaveDM
simulate_waveDM(;
    Nx = Ny = Nz = 64,
    Nt = 3,
    Xmax = 2.0,
    Tmax = 0.001,
    gpu = true,
    Realtime = false,         # headless test
    unicode_plot = true,      # terminal progress
    outputdir = mktempdir(),
)
```

## Enabling high-performance FFT

For multi-threaded FFTW on a CPU node:

```julia
using FFTW
FFTW.set_provider!("mkl")            # or "fftw" if you prefer the bundled backend
FFTW.set_num_threads(Threads.nthreads())
```

Pair this with `julia -t auto` (or `julia -t N`) when launching Julia.

## Enabling distributed memory

Add worker processes from inside your script before loading WaveDM:

```julia
using Distributed
addprocs(4)              # spawn four local workers
@everywhere using WaveDM
```

then pass `distributed_memory = true` to `simulate_waveDM` / `SPE3D_waveDM`. For multi-node clusters use `addprocs([("node1", :auto), ("node2", :auto), ...])` with SSH or a cluster manager.

## Verifying the installation

```julia
pkg> test WaveDM
```

`test WaveDM` runs the unit, integration, and energy-convergence tests in [`test/`](https://github.com/JuliaAstroSim/WaveDM.jl/tree/main/test).

## Recommended workflow

1. Use **Julia ≥ 1.11** (we test against the current LTS and the latest release).
2. **Create a dedicated project environment** for each simulation campaign — reproducibility matters.
3. Launch with **`julia -t auto`** for multi-threading; add `FFTW.set_num_threads(...)` for the FFT path.
4. Use **GPU** for grids $\gtrsim 256^3$; switch back to CPU for very small problems where kernel launch overhead dominates.
5. For multi-node runs, **always profile a small CPU run first** to validate the setup before scaling out.
6. Keep dependencies updated: `pkg> update` regularly.

## Troubleshooting

| Symptom | Likely cause / fix |
| --- | --- |
| `OutOfMemoryError` on GPU | Reduce `Nx`, `Ny`, `Nz` or set `gpu = false` |
| Slow first run | Julia JIT compile time — second runs are much faster; use `PackageCompiler` for a precompiled sysimage |
| `MethodError: no method matching ... ::Float32` | The default `V` in `simulate_waveDM` is `Float32`; supply your own if you need `Float64` |
| `FFTW.PROVIDERS` not picking up MKL | `using FFTW; FFTW.set_provider!("mkl")` must be called **before** the first FFT |
| `Distributed` workers cannot load WaveDM | Add `@everywhere using WaveDM` after `addprocs` |
| `CUDA` initialization fails on Windows | Make sure the CUDA driver matches the CUDA toolkit bundled with your `CUDA.jl` version |
| Profile-fitting stalls | Reduce `extract_min_t` so fewer snapshots are scanned, or pre-fit on a single snapshot |

If you encounter an issue not listed here, please open a [GitHub Issue](https://github.com/JuliaAstroSim/WaveDM.jl/issues) with a minimal reproduction.

## Uninstallation

```julia
pkg> rm WaveDM
```

## Next steps

- Explore the [Examples](@ref) section for end-to-end simulations.
- Read the [Algorithms](@ref) section to understand the SPE, split-step Fourier method, boundary conditions, and QSS diagnostics.
- Check the [API reference](api/configs.md) for the configuration structures.
