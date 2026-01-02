# Utility functions and constants for WaveDM

# Constant parameters from temp/SPE-WaveDM.jl
const kernel3 = ones(3,3,3) ./ 3^3
const kernel5 = ones(5,5,5) ./ 5^3
const kernel7 = ones(7,7,7) ./ 7^3
const kernel9 = ones(9,9,9) ./ 9^3
const kernel11 = ones(11,11,11) ./ 11^3
const kernel15 = ones(15,15,15) ./ 15^3

"""
    filter_min_rho(rho, r, spl)

Check if density `rho` is above the minimum threshold defined by spline `spl` at radius `r`.
"""
function filter_min_rho(rho, r, spl)
    if rho <= spl([r])[1]
        return false
    else
        return true
    end
end

"""
    func_dθ_dt(phase_last_t, phase, dt)

Compute the time derivative of phase using unwrapped phases.
"""
function func_dθ_dt(phase_last_t, phase, dt)
    phase_unwraped = unwrap([phase_last_t, phase])
    return (phase_unwraped[2] - phase_unwraped[1]) / dt
end

"""
    extract_ma_number(s::AbstractString)

Extract the ma number from a string like "ma(3.0)".
"""
function extract_ma_number(s::AbstractString)
    # ma(xx.xx)
    reg = match(r"ma\(([\d.]+)\)", s)
    if reg !== nothing
        return parse(Float64, reg[1])
    else
        return nothing
    end
end
