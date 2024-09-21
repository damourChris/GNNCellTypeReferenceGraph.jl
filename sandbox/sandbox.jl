using Pkg
Pkg.activate(@__DIR__)

using ExpressionData
using EnsemblGeneMapping
using OntologyMapping
using OntologyTrees
using StatsBase

include("../src/GNNCellTypeReferenceGraph.jl")

eset = load_eset(joinpath(@__DIR__, "data", "GSE22886-GPL96_series_matrix.rds"))
mapped_eset = map_to_ensembl(eset, "entrezgene_id"; gene_col="ENTREZ_GENE_ID")

replacements = ("cells" => "cell",
                "Naïve" => "naive",
                "Monocytes" => "monocyte",
                "Macrophages" => "macrophage",
                "Neutrophils" => "neutrophil",
                "CD8+" => "CD8",
                "CD4+" => "CD4",
                "Tregs" => "Treg")

ontology_mapped_eset = map_to_ontology(mapped_eset, "cl"; ontology_col="cell type:ch1",
                                       replacements=replacements)

# Now we need to get the pairings
# -> This step is to establish the links between genes and cell types
# Here, lets just randomly assign some genes to cell types
# while making sure that no gene is assigned to more than one cell type
pairings = Dict{Symbol,Vector{String}}()
for cell_type in unique(phenotype_data(ontology_mapped_eset)[!, "ontology_ids"])
    # generate N random genes
    N = rand(15:50)
    genes = sample(feature_names(ontology_mapped_eset), N; replace=false)
    pairings[Symbol(cell_type)] = genes
end

GNNCellTypeReferenceGraph.iri_from_id("CL:0000000")

# Create the reference graph
onto_tree = GNNCellTypeReferenceGraph.create_reference_graph(pairings)