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
        # the full URL with https://. If not, the host becomes part of the base path
        deploy_url = "https://juliageophysics.com/Timekeepers.jl",
        description = "Timeseries analysis in Electromagnetic Geophysics",
        # dev is the only published version. Thus, let search engines index it
        noindex_non_stable = false,
    ),
    # nested lists become the dropdown menus of the top navigation bar
    pages = [
        "Home"            => "index.md",
        "Getting Started" => "getting_started.md",
        "TKDash"          => "tkdash.md",
        "TKApp"           => "tkapp.md",
        "TKProc"          => "tkproc.md",
        "For Developers" => [
            "Overview"           => "developers/index.md",
            "Instrument Formats" => "developers/formats.md",
            "Masking & Cleaning" => "developers/masking.md",
            "Metronix Sites"     => "developers/metronix.md",
            "TKDash in Detail"   => "developers/tkdash.md",
            "TKApp in Detail"    => "developers/tkapp.md",
            "Spectral Views"     => "developers/spectra.md",
            "API Reference"      => "api.md",
        ],
    ],
    checkdocs = :exports,
    doctest = false,
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/JuliaGeophysics/Timekeepers.jl.git",
    devbranch = "main",
    push_preview = true,
)
