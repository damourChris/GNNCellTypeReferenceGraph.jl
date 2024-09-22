module GNNCellTypeReferenceGraph

include("./setup_env.jl")

using DataFrames
using ExpressionData
using Graphs
using JLD2
using LightXML
using MetaGraphs
using Random
using OntologyLookup
using OntologyTrees
using ProteinMapping
using RCall
using Statistics
using STRINGdb

include("utils.jl")

include("main.jl")
export load_input, create_reference_graph

end
