module GNNCellTypeReferenceGraph

include("./setup_env.jl")

using JLD2
using DataFrames
using Statistics
using Random
using MetaGraphs
using RCall

using ExpressionData
using OntologyLookup
using OntologyTrees

include("utils.jl")

include("main.jl")
export load_input, create_reference_graph

end
