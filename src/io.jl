# SPE input/output module

"""
$(TYPEDSIGNATURES)

Save initial conditions to file.
This function encapsulates the IC saving code from SPE3D_waveDM.
"""
function save_initial_conditions(ψ, baryon_mode, baryon_potential, ax_b, ay_b, az_b, outputdir, title, suffix, particles = nothing)
    if baryon_mode == :particles_dynamic
        save(joinpath(outputdir, "$(title), $(suffix) - IC.jld2"), Dict("ψ" => ψ, "Φ_b" => baryon_potential, "ax_b"=>ax_b, "ay_b"=>ay_b, "az_b"=>az_b, "particles" => particles))
    elseif baryon_mode != :ignored
        save(joinpath(outputdir, "$(title), $(suffix) - IC.jld2"), Dict("ψ" => ψ, "Φ_b" => baryon_potential, "ax_b"=>ax_b, "ay_b"=>ay_b, "az_b"=>az_b))
    else
        save(joinpath(outputdir, "$(title), $(suffix) - IC.jld2"), Dict("ψ" => ψ))
    end
end

"""
$(TYPEDSIGNATURES)

Save simulation evolution results to file.
This function encapsulates the phi saving code from SPE3D_waveDM.
"""
function save_evolution_results(ψ, Φ_all, ψ_last_t, average, averaged_ψ2, outputdir, title, suffix)
    save(joinpath(outputdir, "$(title), $(suffix) - Prop.jld2"),
        "ψ", ψ,
        "Φ_all", Φ_all,
        "ψ_last_t", ψ_last_t,
        "ψ2_averaged", average ? averaged_ψ2 : nothing,
    )
end

"""
$(TYPEDSIGNATURES)

Create and save property DataFrame.
This function encapsulates the DataFrame creation code from SPE3D_waveDM.
"""
function save_property_dataframe(ArrayT, ArrayR, ArrayR1, ArrayR2, ArrayR3, ArrayR4, ArrayR5,
    ArrayR6, ArrayR7, ArrayR8, ArrayR9, ArrayTotalMass, plot_virial,
    ArrayVirialPotential, ArrayTotalKineticE, ArrayTotalQuantumE, ArrayVirial,
    ArrayMomentumX, ArrayMomentumY, ArrayMomentumZ, outputdir, title, suffix)
    
    dfProp = DataFrame(
        :t => ArrayT[],
        :R => ArrayR[],
        :R1 => ArrayR1[],
        :R2 => ArrayR2[],
        :R3 => ArrayR3[],
        :R4 => ArrayR4[],
        :R5 => ArrayR5[],
        :R6 => ArrayR6[],
        :R7 => ArrayR7[],
        :R8 => ArrayR8[],
        :R9 => ArrayR9[],
        :TotalMass => ArrayTotalMass[],
    )
    if plot_virial
        dfProp[!,:PE_abs] = ArrayVirialPotential[]
        dfProp[!,:KE] = ArrayTotalKineticE[]
        dfProp[!,:QE] = ArrayTotalQuantumE[]
        dfProp[!,:Virial] = ArrayVirial[]
        dfProp[!,:MomentumX] = ArrayMomentumX[]
        dfProp[!,:MomentumY] = ArrayMomentumY[]
        dfProp[!,:MomentumZ] = ArrayMomentumZ[]
    end
    CSV.write(joinpath(outputdir, "$(title), $(suffix) - Prop.csv"), dfProp)
    
    return dfProp
end

"""
$(TYPEDSIGNATURES)

Compute averaged density and acceleration fields.
This function encapsulates the averaging code from SPE3D_waveDM.
"""
function compute_averaged_fields(average, buffer_ψ2, average_N, baryon_mode, a_all, Φ_b, Δ, Nx, Ny, Nz, gpu, GPU, CPU, Periodic, fft_poisson, grad_central)
    if average
        averaged_ψ2 = buffer_ψ2 / average_N
        if baryon_mode == :ignored
            averaged_Φ_all = collect(4π * fft_poisson(Δ, [Nx-1, Ny-1, Nz-1], DeviceArray(averaged_ψ2), Periodic(), gpu ? GPU() : CPU()))
        else
            averaged_Φ_all = collect(4π * fft_poisson(Δ, [Nx-1, Ny-1, Nz-1], DeviceArray(averaged_ψ2), Periodic(), gpu ? GPU() : CPU())) + collect(Φ_b)
        end
        averaged_ax_all, averaged_ay_all, averaged_az_all = grad_central(-Δ..., averaged_Φ_all)
        averaged_a_all = sqrt.(averaged_ax_all[:, :, div(end,2)].^2 .+ averaged_ay_all[:, :, div(end,2)].^2 .+ averaged_az_all[:, :, div(end,2)].^2)
    else
        averaged_ψ2 = []
        averaged_a_all = a_all .* 0
    end
    
    return averaged_ψ2, averaged_a_all
end

# Export functions
export save_initial_conditions, save_evolution_results, save_property_dataframe, save_acceleration_dataframe
export compute_averaged_fields
