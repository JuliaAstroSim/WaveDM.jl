# Installation

## Prerequisites

WaveDM.jl requires Julia 1.11 or later. Before installing WaveDM.jl, make sure you have Julia installed on your system. You can download Julia from the [official website](https://julialang.org/downloads/).

### Optional GPU Support

If you plan to use GPU acceleration, you'll need:
- An NVIDIA GPU with CUDA support
- CUDA Toolkit installed (version 11.4 or later recommended)


## Installation Methods

### Using Julia Package Manager

The easiest way to install WaveDM.jl is through the Julia Package Manager. Open Julia and enter the package mode by pressing `]`, then run:

```julia
pkg> add WaveDM
```

### From GitHub Repository

To install the latest development version directly from GitHub:

```julia
pkg> add https://github.com/JuliaAstroSim/WaveDM.jl.git
```

### For Development

If you plan to modify the code, clone the repository and install it in development mode:

```bash
git clone https://github.com/JuliaAstroSim/WaveDM.jl.git
cd WaveDM.jl
julia --project -e "using Pkg; Pkg.develop(Pkg.PackageSpec(path="."))"
```

## Required Dependencies

WaveDM.jl depends on several Julia packages:


## Optional Dependencies

These dependencies are needed for specific features:

| Package | Description | Feature |
|---------|-------------|---------|
| `CUDA` | GPU programming support | GPU acceleration |
| `Distributed` | Parallel computing support | Distributed memory parallelism |

## Setting Up GPU Support

### Installing CUDA.jl

If you have a compatible NVIDIA GPU, you can enable GPU acceleration by installing CUDA.jl:

```julia
pkg> add CUDA
```

### Verifying GPU Support

To verify that GPU support is working correctly, run the following Julia code:

```julia
using CUDA

if CUDA.has_cuda()
    println("CUDA is available")
    println("GPU: ", CUDA.name(CUDA.device()))
    println("Compute capability: ", CUDA.capability(CUDA.device()))
else
    println("CUDA is not available")
end
```

## Verifying Installation

After installing WaveDM.jl, you can verify the installation by running CI tests:
```julia
pkg> test WaveDM
```

## Troubleshooting

### Common Issues

1. **Memory allocation errors**
   - Error: `Out of GPU memory` or similar
   - Solution: Reduce the grid size (decrease `Nx`, `Ny`, `Nz`) or disable GPU acceleration by setting `gpu=false`.
2. **Performance issues**
   - Solution: Ensure you're using a recent Julia version, enable multi-threading with `julia -t auto`, and consider using GPU acceleration for large simulations.

### Getting Help

If you encounter issues not covered here, you can:

1. Check the [GitHub Issues](https://github.com/JuliaAstroSim/WaveDM.jl/issues) page for similar problems
2. Open a new issue on GitHub with a detailed error description and reproduction steps

## Uninstallation

To remove WaveDM.jl from your system:

```julia
pkg> rm WaveDM
```

## Recommended Workflow

For optimal performance and reproducibility, we recommend:

1. Using Julia 1.9 or later
2. Creating a dedicated project environment for each simulation
3. Enabling multi-threading with `julia -t auto`
4. Using GPU acceleration for large simulations (>256^3 grid points)
5. Keeping dependencies updated regularly

## Next Steps

After installing WaveDM.jl, you can:

- Explore the [Examples](@ref Examples) section to learn how to use the package
- Read the [Algorithms](@ref Algorithms) section to understand the underlying methods
- Check the API documentation in the `api` directory for detailed function documentation

Happy simulating!