using Documenter
using HomotopyGetsReal

makedocs(
    sitename = "HomotopyGetsReal.jl",
    modules = [HomotopyGetsReal],
    authors = "Juan Camilo Gonzalez",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://JuanCaGC.github.io/HomotopyGetsReal.jl",
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
    ],
    checkdocs = :exports,
    # Existing docstrings use `[`Foo`](@ref)` for flat-module / private helpers
    # (e.g. `Topology.decompose_1d_curve`, `_track_path_segment!`) that are not
    # Documenter pages. Do not rewrite those docstrings here — warn only.
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(
    repo = "github.com/JuanCaGC/HomotopyGetsReal.jl.git",
    push_preview = true,
)
