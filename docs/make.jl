using Documenter, DocumenterVitepress
using Timekeepers

DocMeta.setdocmeta!(Timekeepers, :DocTestSetup, :(using Timekeepers); recursive = true)

makedocs(;
    modules  = [Timekeepers],
    sitename = "Timekeepers.jl",
    authors  = "JuliaGeophysics community, Pankaj K Mishra, and contributors",
    format   = DocumenterVitepress.MarkdownVitepress(;
        repo       = "github.com/JuliaGeophysics/Timekeepers.jl",
        devbranch  = "main",
        devurl     = "dev",
        # full URL with https://, otherwise the host is taken as part of the base path
        deploy_url = "https://juliageophysics.com/Timekeepers.jl",
        description = "Timeseries analysis in Electromagnetic Geophysics",
        # dev is the only published version, so let search engines index it
        noindex_non_stable = false,
    ),
    # nested lists become the dropdown menus of the top navigation bar
    pages = [
        "Home"            => "index.md",
        "Getting Started" => "getting_started.md",
        "Data" => [
            "Instrument Formats" => "formats.md",
            "Masking & Cleaning" => "masking.md",
            "Metronix Sites"     => "metronix.md",
        ],
        "TKApp" => [
            "TKApp Explorer" => "tkapp.md",
            "Spectral Views" => "spectra.md",
        ],
        "API" => "api.md",
    ],
    checkdocs = :exports,
    doctest = false,
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/JuliaGeophysics/Timekeepers.jl.git",
    devbranch = "main",
    push_preview = true,
)
