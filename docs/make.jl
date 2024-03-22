using NetSurvival
using Documenter

DocMeta.setdocmeta!(NetSurvival, :DocTestSetup, :(using NetSurvival); recursive=true)

makedocs(;
    modules=[NetSurvival],
    authors="Oskar Laverny <oskar.laverny@univ-amu.fr> and contributors",
    sitename="NetSurvival.jl",
    format=Documenter.HTML(;
        canonical="https://JuliaSurv.github.io/NetSurvival.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaSurv/NetSurvival.jl",
    devbranch="main",
)
