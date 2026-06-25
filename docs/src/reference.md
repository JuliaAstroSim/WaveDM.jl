# References

This page collects the scientific literature and software packages cited in the WaveDM.jl paper and documentation. It is meant as a quick reference for users who want to follow up on the theoretical background or compare with alternative implementations.

## 1. Wave Dark Matter Theory

- **Hu, W., Barkana, R. & Gruzberg, A.** (2000). *Fuzzy cold dark matter: the wave properties of ultralight particles.* Phys. Rev. Lett. **85**, 1158. — `hu2000fuzzy`. The original proposal of fuzzy / wave CDM as a way to address small-scale CDM challenges.
- **Hui, L.** (2021). *Wave dark matter.* Annu. Rev. Astron. Astrophys. **59**, 247. — `hui2021wave`. A modern review of wave CDM: de Broglie scale, solitons, granules, observations.
- **Hui, L., Ostriker, J. P., Tremaine, S. & Witten, E.** (2017). *Ultralight scalars as cosmological dark matter.* Phys. Rev. D **95**, 043541. — `hui2017ultralight`. Defines the quantum pressure and "quantum energy" used in virial diagnostics.
- **Chavanis, P.-H.** (2019). *Predictive model of BEC dark matter halos with a solitonic core and a quantum pressure bulge.* Phys. Rev. D **100**, 083022. — `chavanis2019predictive`. Critical-temperature estimate for ultralight bosons.
- **Kagan, Yu., Muryshev, A. E. & Shlyapnikov, G. V.** (1997). *Evolution of correlation properties and appearance of interference in the course of Bose-Einstein condensation.* Phys. Rev. Lett. **81**, 933. — `kagan1997evolutioncorrelation`. Justifies the classical-field GPE in the highly-degenerate regime.
- **Dong, X.-B. et al.** (2025). *The promise of wave dark matter.* (in prep.). — `dong2025promise`. Discussion of the excited-condensate component and its observational consequences (JWST).
- **Liu, W.-J. et al.** (2025). *Warm fuzzy dark matter.* (in prep.). — `liu2025warmfuzzy`. Wave warm dark matter, an alternative to the cold case.
- **Chiang, B. T., Ostriker, J. P. & Schive, H.-Y.** (2021). *Soliton oscillations and revised constraints from the soliton–halo relation and from pulsar timing.* Mon. Not. R. Astron. Soc. **508**, 3945. — `chiang2021soliton`. Ground-state oscillation period $\tau_{00}$.

## 2. Schrödinger–Poisson Equation

- **Pitaevskii, L. P. & Stringari, S.** (2016). *Bose–Einstein Condensation and Superfluidity.* Oxford University Press. — `pitaevskii2016boseeinstein`. Background on the GPE.
- **Liao, S. et al.** (2025). *Deciphering the soliton–halo relation and the quantum pressure in fuzzy dark matter halos.* (in prep.). — `liao2025decipheringsolitonhalo`. Discusses the velocity-dispersion QSS diagnostic and the "effective thermal energy" terminology.
- **Mocz, P. et al.** (2017). *Galaxy formation with BECDM — I. Turbulence and relaxation of idealized halos.* Mon. Not. R. Astron. Soc. **471**, 4559. — `mocz2017galaxy`. Defines the "gradient energy" $K_\rho$ notation used in the virial-diagnostic literature.

## 3. Numerical Methods

- **Edwards, F., Kendall, E., Hotchkiss, S. & Sherwin, R.** (2018). *PyUltraLight: a pseudo-spectral solver for ultralight dark matter dynamics.* (ascl:1808.008). — `edwards2018pyultralight`. The original pseudo-spectral SPE solver in Python; WaveDM.jl follows the same dimensionless-unit convention.
- **Glennon, N., Musoke, N., Chatterjee, A. & Nikoleyko, D.** (2021). *Modifying SPH for simulations of ultralight bosonic dark matter.* Phys. Rev. D **105**, 123540. — `glennon2021modifying`. Mocks the same dimensionless units.
- **Schwabe, B., Gosenca, M., Behrendt, C., Easther, R. & McDonald, J.** (2016). *Simulations of ultralight axion-like dark matter.* Phys. Rev. D **95**, 083513. — `schwabe2016simulations`. CFL timestep criterion.
- **May, S. & Springel, V.** (2021). *Structure formation in large-volume cosmological simulations of fuzzy dark matter.* Mon. Not. R. Astron. Soc. **506**, 2603. — `may2021structure`. Moving-mesh SPE.
- **May, S., Springel, V. & White, S. D. M.** (2023). *The halo structure of fuzzy dark matter.* (in prep.). — `may2023halo`. Velocity resolution criterion.
- **Li, X., Hui, L. & Yavetz, T.** (2019). *Numerical challenge of the fuzzy dark matter model.* Phys. Rev. D **103**, 023508. — `li2019numerical`. Resolution criteria; de Broglie-wavelength requirement.
- **Schive, H.-Y., Chiueh, T. & Broadhurst, T.** (2014, 2026 update). *Cosmic structure as the quantum interference of a coherent dark wave.* Nature Physics **10**, 496; updated review preprint 2026. — `schive2014cosmicstructure`, `schive2026fuzzydark`. The reference paper for the SPE + AMR approach and a comprehensive review of methodologies.
- **James, R. A.** (1977). *The solution of Poisson's equation for isolated source distributions.* J. Comput. Phys. **25**, 71. — `james1977solutionpoissons`. The four-step algorithm used for `boundary = Vacuum()`.
- **Dutta Chowdhury, D., van den Bosch, F. C., Robles, V. H., van Dokkum, P., Schive, H.-Y., Chiueh, T. & Broadhurst, T.** (2021). *On the random motion of the soliton in fuzzy dark matter halos.* Astrophys. J. **916**, 27. — `duttachowdhury2021random`. Box-size and resolution discussion for isolated halos.

## 4. Self-Interaction & Nonlinear Schrödinger

- **Levkov, D., Panin, A. & Tkachev, I.** (2018). *Gravitational Bose-Einstein condensation in the kinetic regime.* Phys. Rev. Lett. **121**, 151301. — `levkov2018gravitational`.
- **Nori, M. & Baldi, M.** (2018). *AX-GADGET: a new code for cosmological simulations of fuzzy dark matter and axion models.* Mon. Not. R. Astron. Soc. **478**, 3935. — `nori2018axgadget`.

## 5. Hydrodynamic / Madelung-based Methods

- **Mocz, P. & Succi, S.** (2015). *Numerical solution of the nonlinear Schrödinger equation using smoothed-particle hydrodynamics.* Phys. Rev. E **91**, 053304. — `mocz2015numerical`.
- **Veltmaat, J., Niemeyer, J. C. & Schwabe, B.** (2016). *Formation and structure of ultralight bosonic dark matter halos.* Phys. Rev. D **94**, 123523. — `veltmaat2016cosmological`.
- **Veltmaat, J., Schwabe, B. & Niemeyer, J. C.** (2018). *Formation of ultralight bosonic dark matter halos through a real scalar field.* Phys. Rev. D **98**, 043509. — `veltmaat2018formation`.
- **Zhang, J., Liu, H. & Chu, M.-C.** (2018). *Axion-Gadget: a smoothed-particle-hydrodynamics code for ultralight dark matter.* (preprint). — `zhang2018ultralight`.
- **Mocz, P. et al.** (2017). *Galaxy formation with BECDM — I. Turbulence and relaxation of idealized halos.* Mon. Not. R. Astron. Soc. **471**, 4559. — `mocz2017galaxy`.
- **Mocz, P. et al.** (2020). *Galaxy formation with BECDM — II. Cosmic filaments and halos.* Mon. Not. R. Astron. Soc. **494**, 2027. — `mocz2020galaxy`.
- **Mocz, P. et al.** (2025). *jaxion: a high-performance code for ultralight dark matter.* (preprint). — `mocz2025jaxion`.
- **Schwabe, B. & Niemeyer, J. C.** (2020). *Simulations of solitonic core plus fuzzy dark halo formation.* Phys. Rev. Lett. **124**, 201301. — `schwabe2020simulating`.
- **Schwabe, B., Gosenca, M., Behrendt, C., Easther, R. & McDonald, J.** (2022). *Deep learning for freckles in BEC dark matter.* Phys. Rev. D **105**, 123502. — `schwabe2022deep`.
- **Lague, A., Schwabe, B., Clough, K. & Niemeyer, J. C.** (2021). *Evolving fuzzy dark matter perturbations.* Phys. Rev. D **104**, 083536. — `lague2021evolving`.
- **Hopkins, P. F. et al.** (2019). *Stable hierarchical solar system and ultra-light dark matter.* Mon. Not. R. Astron. Soc. **488**, 4168. — `hopkins2019stable`.
- **Mina, M., Mota, D. F. & Winther, H. A.** (2020). *SCALAR: an AMR code for fuzzy dark matter.* (preprint). — `mina2020scalar`.
- **Nori, M., Mocz, P., Davé, R. & Baldi, M.** (2024). *Fuzzy-Gasoline: a multi-fluid SPH code for ultralight bosonic dark matter.* (preprint). — `nori2024fuzzy`.

## 6. Related Simulation Codes

- **Musoke, N., Hotchkiss, S. & Easther, R.** (2024). *UltraDark.jl: a Julia package for ultralight dark matter dynamics.* J. Open Source Softw. — `musoke2024ultradarkjl`. Closest Julia-language SPE cousin.
- **Kunkel, S. et al.** (2025). *Hybrid CPU/GPU and MPI parallelization of the GAMER code.* Astrophys. J. Suppl. Ser. — `kunkel2025hybrid`. The MPI/AMR SPE implementation in GAMER.
- **Hlozek, R. et al.** (2015). *A search for ultralight axions using the CMB.* Phys. Rev. D **91**, 103512. — `hlozek2015search`. AxionCAMB.

The full list with method, boundary, capabilities is reproduced in Table 1 of the accompanying paper.

## 7. Observational Data & Galaxy Models

- **Battaglia, G., Taibi, S., Thomas, G. F. & Fritz, T. K.** (2022). *Gaia early DR3 systemic motions of Local Group dwarf galaxies.* Astron. Astrophys. **657**, A89. — `battaglia2022gaia`. Source catalogue for `AstroIC.load_data_MW_satellites()`.
- **Hayashi, K. et al.** (2023). *Dark matter halos of dwarf galaxies from cosmological simulations.* (preprint). — `hayashi2023dark`. Crater II gNFW fit used in the case study.
- **Zhu, Y. et al.** (2023). *How baryonic processes affect the Milky Way mass model.* Mon. Not. R. Astron. Soc. — `zhu2023how`. Multi-component MW mass model (bulge, disk, gas, halo).

## 8. QSS / Quasi-Stationary State Theory

- **Lynden-Bell, D.** (1967). *Statistical mechanics of violent relaxation in stellar systems.* Mon. Not. R. Astron. Soc. **136**, 101. — `lynden-bell1967statisticalmechanics`.
- **Benetti, F. P. C. & Marcos, B.** (2014). *Non-equilibrium phase transitions in a 3D self-gravitating system.* (preprint). — `benetti2014nonequilibrium`.
- **Klessen, R. S. (2000)**. *Gravitational turbulence.* Habilitation thesis. — `klessen2000gravitational`.
- **Berman, R. H. & Tajima, T.** (2012). *Turbulence in self-gravitating systems.* (preprint). — `berman2012turbulence`.
- **Moon, S. et al.** (2024). *Theory of soliton–halo co-evolution in fuzzy dark matter.* (preprint). — `moon2024theory`.
- **Chandrasekhar, S.** (1953). *Problems of stability in gravitational systems.* Rev. Mod. Phys. — `chandrasekhar1953problems`.
- **Woo, T.-P. & Chiueh, T.** (2009). *High-resolution simulation on structure formation with extremely light bosonic dark matter.* Astrophys. J. **697**, 850. — `woo2009highresolution`.

## 9. Software References

### 9.1 Key dependencies

- **Bezanson, J., Edelman, A., Karpinski, S. & Shah, V. B.** (2017). *Julia: a fresh approach to numerical computing.* SIAM Review **59**, 65. — The Julia language.
- **Frigo, M. & Johnson, S. G.** (2005). *The design and implementation of FFTW3.* Proc. IEEE **93**, 216. — `FFTW.jl`.
- **Besard, T., Foket, C. & De Sutter, B.** (2019). *Effective extensible programming: unleashing Julia on GPUs.* IEEE Trans. Parallel Distrib. Syst. — `CUDA.jl`.
- **Byrne, S., Wilcox, L. C. & Hamrud, M.** (2021). *MPI.jl: Julia bindings for the Message Passing Interface.* (preprint). — `MPI.jl`. (Not currently used; planned for a future release.)
- **Danisch, S. & Krumbiegel, J.** (2021). *Makie.jl: flexible high-performance data visualization for Julia.* J. Open Source Softw. — `danisch2021makiejlflexible`. `GLMakie.jl`, `CairoMakie.jl`.
- **Darema, F. et al.** (1988). *A single-program-multiple-data computational model for EPEX/FORTRAN.* Parallel Comput. — `darema1988singleprogrammultipledatacomputational`. The SPMD paradigm underlying WaveDM.jl's `Distributed.jl` integration.
- **Springel, V.** (2005). *The cosmological simulation code GADGET-2.* Mon. Not. R. Astron. Soc. **364**, 1105. — `springel2005cosmological`. Tree / SPH reference.
- **Springel, V., Yoshida, N. & White, S. D. M.** (2001). *GADGET: a simulation environment for cosmological and galactic dynamics.* New Astron. **6**, 79. — `springel2001gadget`.
- **Gropp, W., Lusk, E. & Skjellum, A.** (1999). *Using MPI: portable parallel programming with the message-passing interface.* MIT Press. — `gropp1999usingmpi`.

### 9.2 JuliaAstroSim Ecosystem

- [`PhysicalParticles.jl`](https://github.com/JuliaAstroSim/PhysicalParticles.jl) — vector algebra and N-body particle types.
- [`PhysicalFFT.jl`](https://github.com/JuliaAstroSim/PhysicalFFT.jl) — FFT primitives for physics.
- [`PhysicalFDM.jl`](https://github.com/JuliaAstroSim/PhysicalFDM.jl) — finite-difference operators.
- [`PhysicalMeshes.jl`](https://github.com/JuliaAstroSim/PhysicalMeshes.jl) — mesh bookkeeping.
- [`PhysicalTrees.jl`](https://github.com/JuliaAstroSim/PhysicalTrees.jl) — tree gravity.
- [`AstroSimBase.jl`](https://github.com/JuliaAstroSim/AstroSimBase.jl) — common types and units.
- [`AstroNbodySim.jl`](https://github.com/JuliaAstroSim/AstroNbodySim.jl) — N-body glue and parallel runtime.
- [`AstroIC.jl`](https://github.com/JuliaAstroSim/AstroIC.jl) — initial condition generation.
- [`AstroIO.jl`](https://github.com/JuliaAstroSim/AstroIO.jl) — snapshot I/O.
- [`AstroPlot.jl`](https://github.com/JuliaAstroSim/AstroPlot.jl) — analysis and visualization.
- [`GalacticDynamics.jl`](https://github.com/JuliaAstroSim/GalacticDynamics.jl) — density profiles and orbits.

The full ecosystem is at <https://github.com/JuliaAstroSim>.

## 10. WaveDM.jl Citation

If you use WaveDM.jl in your research, please cite it as:

```bibtex
@article{MengAndDong2026WaveDMjl,
    title = {WaveDM.jl: An Adaptable Simulation Framework for Dynamics of Baryonic and Wave Dark Matter on Galaxy Scales},
   author = {Run-Yu Meng and Xiao-Bo Dong},
  journal = {arXiv preprint},
     year = 2026,
    month = jun,
      url = {https://arxiv.org/abs/2606.25026v1},
  urldate = {2026-06-25}
}
```

## 11. Online Resources

- WaveDM.jl repository — <https://github.com/JuliaAstroSim/WaveDM.jl>
- WaveDM.jl documentation — <https://JuliaAstroSim.github.io/WaveDM.jl/dev>
- JuliaAstroSim organization — <https://github.com/JuliaAstroSim>
- Julia language — <https://julialang.org/>

