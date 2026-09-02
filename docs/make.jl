using Documenter
using Timekeepers

DocMeta.setdocmeta!(Timekeepers, :DocTestSetup, :(using Timekeepers); recursive = true)

makedocs(;
    modules  = [Timekeepers],
    sitename = "Timekeepers.jl",
    authors  = "JuliaGeophysics community, Pankaj K Mishra, and contributors",
    format   = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical  = "https://juliageophysics.github.io/Timekeepers.jl",
        assets     = ["assets/timekeepers.css"],
        collapselevel = 2,
    ),
    pages = [
        "Home"              => "index.md",
        "Getting Started"   => "getting_started.md",
        "Instrument Formats" => "formats.md",
        "Masking & Cleaning" => "masking.md",
        "TKApp Explorer"    => "tkapp.md",
        "Spectral Views"    => "spectra.md",
        "Metronix Sites"    => "metronix.md",
        "API Reference"     => "api.md",
    ],
    checkdocs = :exports,
    doctest = false,
)

deploydocs(;
    repo = "github.com/JuliaGeophysics/Timekeepers.jl.git",
    devbranch = "main",
    push_preview = true,
)
