using GNNCellTypeReferenceGraph
using Documenter

using Changelog

DocMeta.setdocmeta!(GNNCellTypeReferenceGraph, :DocTestSetup,
                    :(using GNNCellTypeReferenceGraph); recursive=true)

Changelog.generate(Changelog.CommonMark(),
                   joinpath(@__DIR__, "../CHANGELOG.md");
                   repo="damourChris/GNNCellTypeReferenceGraph.jl")

Changelog.generate(Changelog.Documenter(),                 # output type
                   joinpath(@__DIR__, "../CHANGELOG.md"),  # input file
                   joinpath(@__DIR__, "src/CHANGELOG.md"); # output file
                   repo="damourChris/GNNCellTypeReferenceGraph.jl",)

makedocs(;
         modules=[GNNCellTypeReferenceGraph],
         authors="Chris Damour",
         sitename="GNNCellTypeReferenceGraph.jl",
         format=Documenter.HTML(;
                                canonical="https://damourChris.github.io/GNNCellTypeReferenceGraph.jl",
                                edit_link="main",
                                assets=String[],),
         pages=["Home" => "index.md",
                "Changelog" => "CHANGELOG.md"],)

deploydocs(;
           repo="github.com/damourChris/GNNCellTypeReferenceGraph.jl",
           devbranch="main",)
