# Filter and utility functions for WaveDM

"""
$(TYPEDSIGNATURES)

Check if density `rho` is above the minimum threshold defined by spline `spl` at radius `r`.

`spl` is expected to be a callable that accepts a vector of radii and returns a vector of threshold values (e.g. a `Dierckx.Spline1D`). 
The function returns `true` when `rho` is strictly greater than the interpolated threshold at `r`.
"""
function filter_min_rho(rho, r, spl)
    if rho <= spl([r])[1]
        return false
    else
        return true
    end
end

"""
$(TYPEDSIGNATURES)

Compute the time derivative of phase using unwrapped phases.
"""
function func_dθ_dt(phase_last_t, phase, dt)
    phase_unwraped = unwrap([phase_last_t, phase])
    return (phase_unwraped[2] - phase_unwraped[1]) / dt
end

"""
$(TYPEDSIGNATURES)

Extract the ma number from a string like "ma(3.0)".
"""
function extract_ma_number(s::AbstractString)
    reg = match(r"ma\(([\d.]+)\)", s)
    if reg !== nothing
        return parse(Float64, reg[1])
    else
        return nothing
    end
end
