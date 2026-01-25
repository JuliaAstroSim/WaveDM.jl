import Base.length
import Base.iterate
@inline length(spl::Spline1D) = 1
@inline iterate(spl::Spline1D) = (spl,nothing)
@inline iterate(spl::Spline1D,st) = nothing

"""
$(TYPEDSIGNATURES)

Interpolate data `(x, v)` at points in `xp`

Use `Dierckx.Spline1D`:
- Default k = 1 linear interpolation might be better in signal processing
- k = 2 quadratic
- k = 3 cubic
- k <= 5

!!! attentions
    - Data must be unitless, and have no NaN
    - Coordinate array is arranged from small to large
"""
function interp1(x, v, xp; k = 1, kw...)
    spl = Spline1D(ustrip.(unit(eltype(x)), x), v; k, kw...)
    return spl(ustrip.(unit(eltype(xp)), xp))
end

"""
$(TYPEDSIGNATURES)

Mean Squared Error (MSE)
"""
function mse(y_pred, y_true)
    return mean((y_pred .- y_true) .^ 2)
end

chi2(obs, sigma, pred) = sum(((obs.-pred)./sigma).^2)

chi2reduced(obs, sigma, pred; DOF=0) = sum(((obs.-pred)./sigma).^2) / (length(pred)-DOF)

chi2reduced(obs, sigma, pred, err; DOF=0) = sum((obs.-pred).^2 ./ (sigma.^2 + err.^2)) / (length(pred)-DOF)

"""
$(TYPEDSIGNATURES)

χ²

The observation and prediction are not one-to-one.
Interpolate the discrete observational data using 1D spline and compute χ² within observational region
"""
function chi2(obs, sigma, obs_r, pred, r; 
    k = 1, # 1, linear; 2, quadratic; 3, cubic
    kw...
)
    rMin = minimum(obs_r)
    rMax = maximum(obs_r)
    indices = rMin .<= r .<= rMax
    return chi2(
        interp1(obs_r, obs, r[indices]; k, kw...),
        interp1(obs_r, sigma, r[indices]; k, kw...),
        pred[indices],
    )
end

"""
$(TYPEDSIGNATURES)
reduced χ²

The observation and prediction are not one-to-one.
Interpolate the discrete observational data using 1D spline and compute χ² within observational region
"""
function chi2reduced(obs, sigma, obs_r, pred, r; 
    k = 1, # 1, linear; 2, quadratic; 3, cubic
    DOF = 0,
    kw...
)
    rMin = minimum(obs_r)
    rMax = maximum(obs_r)
    indices = rMin .<= r .<= rMax
    return chi2reduced(
        interp1(obs_r, obs, r[indices]; k, kw...),
        interp1(obs_r, sigma, r[indices]; k, kw...),
        pred[indices];
        DOF,
    )
end

"""
$(TYPEDSIGNATURES)
reduced χ² with erros from prediction

The observation and prediction are not one-to-one.
Interpolate the discrete observational data using 1D spline and compute χ² within observational region
"""
function chi2reduced(obs, sigma, obs_r, pred, err, r; 
    k = 1, # 1, linear; 2, quadratic; 3, cubic
    DOF = 0,
    kw...
)
    rMin = minimum(obs_r)
    rMax = maximum(obs_r)
    indices = rMin .<= r .<= rMax
    return chi2reduced(
        interp1(obs_r, obs, r[indices]; k, kw...),
        interp1(obs_r, sigma, r[indices]; k, kw...),
        pred[indices],
        err[indices];
        DOF,
    )
end


function velocity_anisotropy(vx, vy, vz)
    σ = [std(vx), std(vy), std(vz)]
    max_σ = maximum(σ)
    min_σ = minimum(σ)
    return anisotropy = 1 - (min_σ / max_σ) # 0 for isotropy, 1 for anisotropy
end


function find_first_intersection(x, y, p)
    for i in 1:length(x)-1
        x1, x2 = x[i], x[i+1]
        y1, y2 = y[i], y[i+1]
        
        # compute intersection point
        if (y1 - p) * (y2 - p) <= 0
            t = (p - y1) / (y2 - y1)  # ratio
            intersection = x1 + t * (x2 - x1)
            return intersection
        end
    end
    
    return nothing # no intersection
end

function find_intersections(x, y, p)
    intersections = Float64[]
    
    for i in 1:length(x)-1
        x1, x2 = x[i], x[i+1]
        y1, y2 = y[i], y[i+1]
        
        # compute intersection point
        if (y1 - p) * (y2 - p) <= 0
            t = (p - y1) / (y2 - y1)  # ratio
            intersection = x1 + t * (x2 - x1)
            push!(intersections, intersection)
        end
    end
    
    return intersections
end
