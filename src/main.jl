
function load_input(eset_file_path::String)

    # We can use the R script to get the cell type marker genes
    # This is a temporary solution until we have a better way to get the cell type marker genes
    r_script_path = joinpath(@__DIR__, "get_cell_type_marker_genes.R")
    eset = load_eset(eset_file_path)

    @rput eset r_script_path
    R"""
    source(r_script_path)

    pairings <- get_cell_type_marker_genes(eset)
    """

    pairs = @rget pairings

    # Convert the pairings to a Julia dictionary
    ref_eset_pairings = Dict{Symbol,Vector{String}}()
    for (key, value) in pairings_R
        ref_eset_pairings[key] = value
    end

    return ref_eset_pairings
end

function create_reference_graph(pairings::Dict{Symbol,Vector{String}}; onto="cl",
                                base_term_iri="http://purl.obolibrary.org/obo/CL_0000000")::OntologyTree
    base_term = onto_term(onto, base_term_iri)
    required_terms_ids = collect(keys(pairings))

    # Required nodes for the graph
    # Here for now we are going to use the first nodes to be Terms (specificaly CL terms)
    term_paired_nodes = Dict{Term,Vector{String}}()
    required_terms = Term[]
    for term_id in sort(required_terms_ids)
        term = onto_term(onto, iri_from_id(term_id))
        if ismissing(term)
            @warn "Term not found: $term_id"
            continue
        end
        push!(required_terms, term)
        # Add the term to the pairings object
        term_paired_nodes[term] = pairings[term_id]
    end

    onto_tree = OntologyTree(base_term, required_terms; max_parent_limit=20)

    add_genes!(onto_tree, unique(vcat(values(term_paired_nodes)...)))

    connect_term_genes!(onto_tree, term_paired_nodes)

    return onto_tree
end