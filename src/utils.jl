function iri_from_id(id::String)::String
    base_url = "http://purl.obolibrary.org/obo/"
    # If the id is in the form with the semi-colon
    # Then we need to replace it with the underscore to get the IRI url 
    if occursin(":", id)
        return "$base_url$(replace(id, ':' => '_'))"
    else
        return "$base_url$id"
    end
end
function iri_from_id(id::Symbol)::String
    return iri_from_id(string(id))
end
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

function split_array(arr, L)
    K = length(arr)

    if K < L
        return [arr]
    end

    n = div(K - 1, L) + 1
    r = mod(K, L)

    result = Vector{Vector{eltype(arr)}}()

    start_idx = 1
    for i in 1:(n - 1)
        push!(result, view(arr, start_idx:(start_idx + L - 1)))
        start_idx += L
    end

    if r > 0
        push!(result, view(arr, start_idx:(start_idx + r - 1)))
    end

    return result
end

function get_pep_mapping(ids::Vector{String})
    df = DataFrame()

    for id_chunk in split_array(ids, 49)
        result = map_to_stable_ensembl_peptide(id_chunk, "ensembl_gene_id")
        append!(df, result)
    end

    return df
end

function get_interactions_batched(ids::Vector{String})
    df = DataFrame()
    chunks = split_array(ids, 50)

    for (index, id_chunk) in enumerate(chunks)
        println("Processing chunk $index / $(length(chunks))")
        result = STRINGdb.get_interactions(id_chunk)
        append!(df, result)
    end

    return df
end

function export_to_graphxml(tree::OntologyTree, filename::String)
    graph = tree.graph
    xdoc = XMLDocument()

    # Create the root element
    xroot = create_root(xdoc, "graphml")
    set_attribute(xroot, "xmlns", "http://graphml.graphdrawing.org/xmlns")

    # For each fieldnames in the Term struct, create a key element
    for field in fieldnames(Term)
        xkey = new_child(xroot, "key")
        set_attribute(xkey, "id", field)
        set_attribute(xkey, "for", "node")
        set_attribute(xkey, "attr.name", field)
        set_attribute(xkey, "attr.type", "string")
    end

    for field in [:id, :type, :label]
        xkey = new_child(xroot, "key")
        set_attribute(xkey, "id", field)
        set_attribute(xkey, "for", "node")
        set_attribute(xkey, "attr.name", field)
        set_attribute(xkey, "attr.type", "string")
    end

    for field in [:expression, :proportion]
        xkey = new_child(xroot, "key")
        set_attribute(xkey, "id", field)
        set_attribute(xkey, "for", "node")
        set_attribute(xkey, "attr.name", field)
        set_attribute(xkey, "attr.type", "float")
    end

    for field in [:weight]
        xkey = new_child(xroot, "key")
        set_attribute(xkey, "id", field)
        set_attribute(xkey, "for", "edge")
        set_attribute(xkey, "attr.name", field)
        set_attribute(xkey, "attr.type", "float")
    end

    # Define graph attributes
    xgraph = new_child(xroot, "graph")
    set_attribute(xgraph, "id", "G")
    set_attribute(xgraph, "edgedefault", "directed") # Assuming directed graph

    vertice_with_term_prop = [v_index
                              for (v_index, v_props) in graph.vprops
                              if v_props[:type] == :term]

    # Add term node metadata
    for v in vertice_with_term_prop
        xnode = new_child(xgraph, "node")
        set_attribute(xnode, "id", string(v))

        term = get_prop(graph, v, :term)

        xdata = new_child(xnode, "data")
        set_attribute(xdata, "key", :type)
        add_text(xdata, "cell")

        for field in fieldnames(Term)
            xdata = new_child(xnode, "data")
            set_attribute(xdata, "key", field)
            add_text(xdata, string(getfield(term, field)))
        end

        try
            proportion = get_prop(graph, v, :proportion)
            if ismissing(proportion)
                continue
            end
            xdata = new_child(xnode, "data")
            set_attribute(xdata, "key", :proportion)
            add_text(xdata, string(proportion))
        catch
            continue
        end
    end

    vertice_with_gene_prop = [v_index
                              for (v_index, v_props) in graph.vprops
                              if v_props[:type] == :gene]

    # Add term node metadata
    for v in vertice_with_gene_prop
        xnode = new_child(xgraph, "node")
        set_attribute(xnode, "id", string(v))

        id = get_prop(graph, v, :id)
        xdata = new_child(xnode, "data")
        set_attribute(xdata, "key", :id)
        add_text(xdata, string(id))

        gene_expression = get_prop(graph, v, :expression)
        if (!ismissing(gene_expression))
            xdata = new_child(xnode, "data")
            set_attribute(xdata, "key", :expression)
            add_text(xdata, string(gene_expression))
        end

        xdata = new_child(xnode, "data")
        set_attribute(xdata, "key", :label)
        add_text(xdata, string(id))

        xdata = new_child(xnode, "data")
        set_attribute(xdata, "key", :type)
        add_text(xdata, "gene")
    end

    # Add edges
    for e in edges(graph)
        xedge = new_child(xgraph, "edge")
        set_attribute(xedge, "source", string(src(e)))
        set_attribute(xedge, "target", string(dst(e)))

        if !has_prop(graph, e, :weight)
            continue
        end

        weight = get_prop(graph, e, :weight)
        xdata = new_child(xedge, "data")
        set_attribute(xdata, "key", "weight")
        add_text(xdata, string(weight))
    end

    # Save XML document to file
    return save_file(xdoc, filename)
end