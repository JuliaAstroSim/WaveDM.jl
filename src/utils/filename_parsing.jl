"""
$(TYPEDSIGNATURES)

Extract `Nx`, `Nt`, `ma` and `Xmax` from a WaveDM output filename.

The expected filename convention is

```
..._Nx(<n>)_..._Nt(<n>)_..._ma(<x>)_..._Xmax(<x>)...jld2
```

Any tag that is not present in the filename is returned as `nothing`.

# Example
```julia
julia> parse_waveDM_filename("run_ma(1.50)_Nx(256)_Nt(1000)_Xmax(20.0).jld2")
(Nx = 256, Nt = 1000, ma = 1.5, Xmax = 20.0)
```
"""
function parse_waveDM_filename(filename::AbstractString)
    reg_Nx   = match(r"Nx\(([\d.]+)\)",   filename)
    reg_Nt   = match(r"Nt\(([\d.]+)\)",   filename)
    reg_ma   = match(r"ma\(([\d.]+)\)",   filename)
    reg_Xmax = match(r"Xmax\(([\d.]+)\)", filename)

    Nx   = reg_Nx   === nothing ? nothing : parse(Int,     reg_Nx[1])
    Nt   = reg_Nt   === nothing ? nothing : parse(Int,     reg_Nt[1])
    ma   = reg_ma   === nothing ? nothing : parse(Float64, reg_ma[1])
    Xmax = reg_Xmax === nothing ? nothing : parse(Float64, reg_Xmax[1])

    return (; Nx, Nt, ma, Xmax)
end
