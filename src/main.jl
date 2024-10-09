function connect_genes!(onto_tree::OntologyTree, pairings::Dict{Symbol,Vector{String}})
    # Get the genes that are in the ontology tree
    graph = onto_tree.graph
    genes_to_map = unique(reduce(vcat, values(pairings)))

    # We assume that the genes are already mapped to Ensembl IDS
    # If not exit 
    # TODO: Conditionally load EnsemblGeneMapping for this
    if !all([startswith(gene, "ENSG") for gene in genes_to_map])
        error("Genes are not mapped to Ensembl IDs")
    end

    @info "Mapping $(length(genes_to_map)) genes to peptides"
    gene_pep_mapping_result = get_pep_mapping(genes_to_map)

    # We can convert this to a dictionary for faster lookup
    # This will have unique peptides as keys and a vector of genes as values
    pep_gene_mapping = Dict{String,Vector{String}}()
    for (gene, pep) in
        zip(gene_pep_mapping_result[!, "ensembl_gene_id"],
            gene_pep_mapping_result[!, "ensembl_peptide_id"])
        if haskey(pep_gene_mapping, pep)
            push!(pep_gene_mapping[pep], gene)
        else
            pep_gene_mapping[pep] = [gene]
        end
    end

    @info "Getting interactions for $(length(keys(pep_gene_mapping))) peptides"
    interactions_result::DataFrame = get_interactions_batched(collect(keys(pep_gene_mapping)))

    # We can convert this to a dictionary for faster lookup
    # Keys: Tuple of peptides
    # Values: Interaction score
    interactions = Dict{Tuple{String,String},Float64}()
    for interaction in eachrow(unique(interactions_result))
        (; stringId_A, stringId_B, score) = interaction

        # The strings will have the format "ORGANISM_ID.PEP_ID"
        # We only want the PEP_ID
        pep1 = string(split(stringId_A, ".")[2])
        pep2 = string(split(stringId_B, ".")[2])

        interactions[(pep1, pep2)] = score
    end

    # Now that have interactions scores, we need to translate back into interactions strength 
    # between genes
    # Keys: Tuple of genes
    # Values: Interaction score
    # Need to make sure that we deal with multiple mappings 
    # -> We can either take the average, max, min, sum or median of the scores
    # -> Here, we will take the average
    # We can also add a threshold to filter out weak interactions
    # For performance, we can check ahead of time what gene-gene interaction have multiple scores and reduce them 
    # Then we can combine everything
    gene_interactions = Dict{Tuple{String,String},Float64}()
    for (pep1, pep2) in keys(interactions)
        # If the interactions has mappings to genes that are not in the gene list, we can skip them
        if !haskey(pep_gene_mapping, pep1) || !haskey(pep_gene_mapping, pep2)
            continue
        end

        genes1 = pep_gene_mapping[pep1]
        genes2 = pep_gene_mapping[pep2]

        for gene1 in genes1, gene2 in genes2
            # Altought we are not considering self-interactions, we should check for them
            # -> This is because the same gene can have multiple peptides and thus multiple interactions
            if gene1 == gene2
                continue
            end

            key = (gene1, gene2)
            if haskey(gene_interactions, key)
                gene_interactions[key] = (gene_interactions[key] +
                                          interactions[(pep1, pep2)]) /
                                         2
            else
                gene_interactions[key] = interactions[(pep1, pep2)]
            end
        end
    end

    # We can now add the edges with weights to the graph
    for (gene1, gene2) in keys(gene_interactions)
        # Find the indinces for both genes
        gene1_idx = graph[gene1, :id]
        gene2_idx = graph[gene2, :id]
        add_edge!(graph, gene1_idx, gene2_idx)

        # Add the weight to the edge
        set_prop!(graph, gene1_idx, gene2_idx, :weight, gene_interactions[(gene1, gene2)])
    end
end

function create_reference_graph(pairings::Dict{Symbol,Vector{String}}; onto="cl",
                                base_term_iri="http://purl.obolibrary.org/obo/CL_0000000",
                                terms_to_exclude=Term[],
                                allow_multiple_root_terms=false,
                                connect_genes::Bool=false)::OntologyTree
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

    onto_tree = OntologyTree(base_term, required_terms, terms_to_exclude;
                             max_parent_limit=20,
                             allow_multiple_roots=false)

    add_genes!(onto_tree, unique(vcat(values(term_paired_nodes)...)))

    if connect_genes
        connect_genes!(onto_tree, pairings)
    end

    connect_term_genes!(onto_tree, term_paired_nodes)

    return onto_tree
end