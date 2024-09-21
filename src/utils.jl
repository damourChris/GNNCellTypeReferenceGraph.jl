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