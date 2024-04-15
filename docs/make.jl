using NetSurvival
using Documenter
using DocumenterCitations

DocMeta.setdocmeta!(NetSurvival, :DocTestSetup, :(using NetSurvival); recursive=true)

bib = CitationBibliography(
    joinpath(@__DIR__,"src","assets","references.bib"),
    style=:numeric
)

makedocs(;
    modules=[NetSurvival],
    plugins=[bib],
    authors="Oskar Laverny <oskar.laverny@univ-amu.fr> and contributors",
    sitename="NetSurvival.jl",
    format=Documenter.HTML(;
        canonical="https://JuliaSurv.github.io/NetSurvival.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "theory.md",
        "getting_started.md",
        "examples.md",
        "references.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaSurv/NetSurvival.jl",
    devbranch="main",
)
