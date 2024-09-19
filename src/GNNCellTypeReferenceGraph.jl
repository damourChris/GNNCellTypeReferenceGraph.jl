module GNNCellTypeReferenceGraph
using JLD2
using DataFrames
using Statistics
using Random
using MetaGraphs
using RCall

using ExpressionData
using OntologyLookup

function __init__(precompile=false)
    return include(joinpath(@__DIR__, "./env.jl"))
end
include("main.jl")
export load_input, create_reference_graph

end
