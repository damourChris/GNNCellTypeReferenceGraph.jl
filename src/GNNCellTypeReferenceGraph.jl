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
    return include("./env.jl")
end

end
