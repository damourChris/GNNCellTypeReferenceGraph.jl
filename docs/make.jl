using GNNCellTypeReferenceGraph
using Documenter

DocMeta.setdocmeta!(GNNCellTypeReferenceGraph, :DocTestSetup, :(using GNNCellTypeReferenceGraph); recursive=true)

makedocs(;
    modules=[GNNCellTypeReferenceGraph],
    authors="Chris Damour",
    sitename="GNNCellTypeReferenceGraph.jl",
    format=Documenter.HTML(;
        canonical="https://damourChris.github.io/GNNCellTypeReferenceGraph.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/damourChris/GNNCellTypeReferenceGraph.jl",
    devbranch="main",
)
