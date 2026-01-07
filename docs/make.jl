"""
Compile with:
julia --project=docs/ --color=yes docs/make.jl

Generate key:
DocumenterTools.genkeys(user="JuliaAstroSim", repo="git@github.com:JuliaAstroSim/PhysicalParticles.jl.git")
"""

using Documenter

using WaveDM

# The DOCSARGS environment variable can be used to pass additional arguments to make.jl.
# This is useful on CI, if you need to change the behavior of the build slightly but you
# can not change the .travis.yml or make.jl scripts any more (e.g. for a tag build).
if haskey(ENV, "DOCSARGS")
    for arg in split(ENV["DOCSARGS"])
        (arg in ARGS) || push!(ARGS, arg)
    end
end


makedocs(
    modules = [WaveDM],
    sitename = "WaveDM.jl",
    authors = "islent",
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://JuliaAstroSim.github.io/WaveDM.jl/dev",
        # assets = ["assets/style.css"],
        analytics = "UA-153693590-1",
        highlights = ["llvm", "yaml"],
    ),
    pages = [
        "Home" => "index.md",
        "Introdution" => "introduction.md",
        "Installation" => "installation.md",
        "Algorithms" => "algorithms.md",
        "APIs" => [
            "Configs" => "api/configs.md",
            # "Solvers" => "api/solver.md",
            # "Initial conditions" => "api/ic.md",
            # "Visualization" => "api/plot.md",
            "KDK" => "api/KDK.md",
            "Simulations" => "api/simulation.md",
            # "APIs" => "api/index.md",
        ],
        "Examples" => "examples.md",
        "References and citations" => "reference.md",
    ],
    clean = false,
    # doctest = true,
    linkcheck = false,
    # linkcheck = true,
    warnonly = [:missing_docs, :cross_references, :docs_block],
)

# deploydocs(
#     repo = "github.com/JuliaAstroSim/WaveDM.jl.git",
#     target = "build",
# )
