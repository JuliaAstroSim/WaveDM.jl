"""
Return Rd in kpc. RHI in kpc
"""
function Rd_from_RHI(RHI, MHI)
    f(x) = log(1e3 * MHI / (2π)) - 2 * log(x) - RHI/x
    
    x = collect(0.01:0.01:50) * RHI
    index = findfirstzero(f.(x))
    if isnothing(index)
        return 0
    else
        return x[index]
    end
end
#TODO test

"""
$(TYPEDSIGNATURES)
"""
function load_SPARC_RC()
    df_SPARC_RC = DataFrame(CSV.File(joinpath(@__DIR__, "../data/SPARC_LTGs/MassModels_Lelli2016c.mrt.txt"), skipto=26, delim=" ", ignorerepeated=true, header=false));
    rename!(df_SPARC_RC, ["Galaxy", "D", "R", "Vobs", "e_Vobs", "Vgas", "Vdisk", "Vbul", "SBdisk", "SBbul"])
    return df_SPARC_RC
end

"""
$(TYPEDSIGNATURES)
"""
function load_SPARC_data(;
    exclude_list = ["F561-1"], #TODO
    prior = "LCDM",
    model = "NFW",
)
    ### Load data
    df_SPARC = DataFrame(CSV.File(joinpath(@__DIR__, "../data/SPARC_LTGs/SPARC_Lelli2016c.mrt.txt"), skipto=99, delim=" ", ignorerepeated=true, header=false))
    rename!(df_SPARC, ["Galaxy", "T", "D", "e_D", "f_D", "Inc", "e_Inc", "L", "e_L", "Reff", "SBeff", "Rdisk", "SBdisk", "MHI", "RHI", "Vflat", "e_Vflat", "Q", "Ref"])

    df_SPARC_bulges = DataFrame(CSV.File(joinpath(@__DIR__, "../data/SPARC_LTGs/Bulges.mrt.txt"), skipto=8, delim="\t", ignorerepeated=true, header=false))
    rename!(df_SPARC_bulges, ["Galaxy", "Lbul"])

    # Halo
    df_SPARC_halo = DataFrame(CSV.File(joinpath(@__DIR__, "../data/SPARC_LTGs/Halo/parameter_$(model)_$(prior).mrt"), skipto=33, delim=" ", ignorerepeated=true, header=false))
    rename!(df_SPARC_halo, ["Galaxy", "Ydisk", "e_Ydisk", "Ybul", "e_Ybul", "D", "e_D", "inc", "e_inc", "V200", "e_V200", "C200", "e_C200", "rs", "e_rs", "log_rhos", "e_log_rhos", "log_M200", "e_log_M200", "Chi"])

    ### Clean data
    # Join by name
    df = innerjoin(select(df_SPARC, Not([:D, :e_D, :f_D, :Ref])), df_SPARC_bulges, select(df_SPARC_halo, Not([:Ybul, :e_Ybul, :C200, :e_C200])), on = :Galaxy)
    
    # Filter out galaxies with asymmetric rotation curves that do not trace the equilibrium gravitational potential (Q = 3)
    filter!(:Q => x->x!=3, df)

    # Filter out those with bulges
    filter!(:Lbul => iszero, df)

    # Filter out those without HI gas
    filter!(:RHI => !iszero, df)

    # Filter out those with zero flat velocity
    filter!(:Vflat => !iszero, df)

    # Filter out face-on galaxies with i < 30◦
    filter!(:Inc => x->x>=30, df)

    # Filter out low luminosity galaxies. Actually filter out UGCA444, which is causing simulation error
    filter!(:L => x->x>=0.1, df)

    # Filter out those with too large halos
    halo_size_max_ratio = 10
    halo_size_min_ratio = 0.1
    filter!([:rs, :RHI] => (rs, RHI) -> halo_size_min_ratio < rs/RHI < halo_size_max_ratio, df)

    #? Exclude galaxies with large error in rotation curves

    # find the radius of Vflat
    df_SPARC_RC = load_SPARC_RC()

    df[!, :Rflat] = zeros(size(df,1))
    for i in 1:size(df,1)
        df_galaxy = filter(:Galaxy => x->x==df.Galaxy[i], df_SPARC_RC)
        index = findvalue(df_galaxy.Vobs, df.Vflat[i])
        if !isnothing(index)
            df.Rflat[i] = df_galaxy.R[first(index)]
        end
    end
    filter!(:Rflat => !iszero, df)

    # solve the scale radius of HI disc
    df[!, :Rd] = Rd_from_RHI.(df.RHI, df.MHI)
    filter!(:Rd => !iszero, df)

    df[!, :a0_mean] = zeros(size(df,1))
    df[!, :a0_std] = zeros(size(df,1))
    df[!, :b_r_IC] = zeros(size(df,1))
    df[!, :b_r] = zeros(size(df,1))
    df[!, :M_b] = zeros(size(df,1))
    df[!, :MHI_sim] = zeros(size(df,1))
    df[!, :M_h] = zeros(size(df,1))
    return df
end

function load_li2018_SPARC()
    df_li2018 = DataFrame(CSV.File(joinpath(@__DIR__, "../data/SPARC_LTGs/Li2018_SPARC_RAR.txt"), delim="\t", ignorerepeated=true, header=2))
    df_li2018.LMR_disk = parse.(Measurements.Measurement{Float64}, df_li2018.LMR_disk)
    df_li2018.LMR_bulge = parse.(Measurements.Measurement{Float64}, df_li2018.LMR_bulge)
    df_li2018.D = parse.(Measurements.Measurement{Float64}, df_li2018.D)
    df_li2018.Inc = parse.(Measurements.Measurement{Float64}, df_li2018.Inc)
    return df_li2018
end