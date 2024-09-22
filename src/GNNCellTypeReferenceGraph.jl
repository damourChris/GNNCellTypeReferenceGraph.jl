module GNNCellTypeReferenceGraph

include("./setup_env.jl")

using JLD2
using DataFrames
using Statistics
using Random
using MetaGraphs
using RCall
using Graphs
using ExpressionData
using OntologyLookup
using OntologyTrees
using ProteinMapping
using LightXML
using STRINGdb

include("utils.jl")

include("main.jl")
export load_input, create_reference_graph

end
