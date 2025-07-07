### -*- Mode: Julia -*-

### algorithms.jl

"""
Find all approximate matches of a query in a DNA sequence.

# Args:
   - `query`: The search query.
   - `tolerance`: Allowed mismatch percentage.
   - `seq`: The LongDNA{4} sequence to search in.

# Returns:
   - A vector of UnitRanges representing match positions.
"""
function approximate_findall(query, tolerance :: Int64, seq :: LongDNA{4})
    results = UnitRange{Int64}[]
    pos = findfirst(query, tolerance, seq)
    while pos !== nothing
        push!(results, pos)
        pos = findnext(query, tolerance, seq, last(pos)+tolerance+1)
    end
    return results
end


"""
Perform approximate search in a JournaledString.

# Args: 
   - `jss`: The JournaledString to search in. 
   - `needle`: The LongDNA{4} sequence to find.

# Returns: 
   - A dictionary mapping time points to match ranges.
"""
function approximate_search(jss :: JournaledString, needle :: LongDNA{4})
    query = ApproximateSearchQuery(needle)
    tolerance = ceil(Int64, (length(needle) / 100) * 5)  
    indexMatrix = Dict{Int64, Vector{UnitRange{Int64}}}()
    vector = approximate_findall(query, tolerance, jss.reference)
    to_remove = Set{UnitRange{Int64}}()
    to_add = Set{UnitRange{Int64}}()

    for i in 1:length(jss.deltaMap)
        indexMatrix[i] = vector
    end
    
    for i in 1:length(jss.deltaMap)
        empty!(to_add)
        empty!(to_remove)
        for range in indexMatrix[i]
            for ( _, entry) in jss.deltaMap[i]
                
                if entry.position in range
                    seq = apply_delta(jss.reference, entry)

                    for element in approximate_findall(query, tolerance, seq)
                        push!(to_add, element)
                    end
                    
                    push!(to_remove, range)
                end

            end
        end
        indexMatrix[i] = filter(x -> all(y -> x != y, to_remove),
                                indexMatrix[i])
        append!(indexMatrix[i], to_add)
        indexMatrix[i] = collect(Set(indexMatrix[i]))
    end
    return indexMatrix
end


"""
Perform approximate search in a JournaledString with tolerance.

# Args:
   - `jss`: The JournaledString to search in.
   - `needle`: The LongDNA{4} sequence to find.
   - `tol`: Allowed mismatch percentage (1-99%).

# Returns:
   - A dictionary mapping time points to match ranges.

# Errors:
   - Throws an error if tolerance is ≤0% or ≥100%.
"""
function approximate_search(jss :: JournaledString,
                            needle :: LongDNA{4},
                            tol :: Int64)

    if tol <= 0 || tol >= 100
        error("Tolerance cannot less or 0% or more than 100%")
    end

    tolerance = ceil(Int64, (length(needle) / 100) * tol) 
    query = ApproximateSearchQuery(needle)
    indexMatrix = Dict{Int64, Vector{UnitRange{Int64}}}()
    vector = approximate_findall(query, tolerance, jss.reference)
    to_remove = Set{UnitRange{Int64}}()
    to_add = Set{UnitRange{Int64}}()

    for i in 1:length(jss.deltaMap)
        indexMatrix[i] = vector
    end

    for i in 1:length(jss.deltaMap)
        empty!(to_add)
        empty!(to_remove)
        for range in indexMatrix[i]
            
            for ( time, entry) in jss.deltaMap[i]
                    
                if entry.position in range
                    seq = apply_delta(jss.reference, entry)

                    for element in approximate_findall(query, tolerance, seq)
                        push!(to_add, element)
                    end
                        
                    push!(to_remove, range)
                end

            end

        end 
        indexMatrix[i] = filter(x -> all(y -> x != y, to_remove),
                                indexMatrix[i])
        append!(indexMatrix[i], to_add)
        indexMatrix[i] = collect(Set(indexMatrix[i]))
    end
    return indexMatrix
end


"""
Perform approximate search in a Journaled String Tree (JST).

# Args:
   - `jst`: The JSTree to search in.
   - `needle`: The LongDNA{4} sequence to find.

# Returns:
    A dictionary mapping node names to match ranges.
"""
function approximate_search(jst :: JSTree, needle :: LongDNA{4})

    seq = jst.root
    query = ApproximateSearchQuery(needle)
    vector = UnitRange{Int}[]
    indexMatrix = Dict{String,
                       Vector{UnitRange{Int64}}}(
                           name => UnitRange{Int64}[] for name in keys(jst.children))
    tolerance = ceil(Int64, (length(needle) / 100) * 5)

    vector = approximate_findall(query, tolerance, seq)
    to_remove = Set{UnitRange{Int64}}()
    to_add = Set{UnitRange{Int64}}()
    range2 = Set{UnitRange{Int64}}()

    for (name, _) in indexMatrix
        indexMatrix[name] = vector
    end

    for (name, node) in jst.children

        for range in indexMatrix[name]

            empty!(to_add)
            empty!(to_remove)
            if !isnothing(node.deltaMap)
                empty!(range2)
                seq = flatten(jst, node.parent.name)
                range2 = push!(range2, range)
                range2 = union!(range2, approximate_findall(query, 
                tolerance, seq))
                union!(to_add, range2)
                for (time, entry) in node.deltaMap
                    for range3 in range2
                        if entry.position in range3
                            
                            seq = apply_delta(seq, entry)
                                
                            for element in
                                approximate_findall(query, tolerance, seq)
                                push!(to_add, element)
                            end
                            push!(to_remove, range3)   
                        end
                    end
                end
                append!(indexMatrix[name], to_add)
                indexMatrix[name]= filter(x -> all(y -> x != y, to_remove), 
                                          indexMatrix[name])
                
            end
            indexMatrix[name] = collect(Set(indexMatrix[name]))
        end
    end
    return indexMatrix
end


"""
Perform exact search in a JournaledString.

# Args:
   - `jss`: The JournaledString to search in.
   - `needle`: The LongDNA sequence to find.

# Returns:
   - A dictionary mapping time points to exact match ranges.
"""
function exact_search(jss :: JournaledString, needle :: LongDNA )
    results = Dict(i => UnitRange{Int64}[] for i in 1:length(jss.deltaMap))
    query = ExactSearchQuery(needle)
    vector = UnitRange{Int64}[]
        for i in 1:length(jss.deltaMap)

            seq = apply_delta(jss.reference, jss.deltaMap[i])
            vector = BioSequences.findall(query, seq)
            append!(results[i], vector)

        end
    return results
end


"""
Perform exact search in a Journaled String Tree (JST).

# Args: 
   - `jst`: The JSTree to search in. 
   - `needle`: The LongDNA{4} sequence to find.

# Returns: 
   - A dictionary mapping node names to exact match ranges.
"""
function exact_search(jst :: JSTree, needle :: LongDNA{4})
    indexMatrix = Dict{String, Vector{UnitRange{Int64}}}()
    query = ExactSearchQuery(needle)
    vector = UnitRange{Int}[]

    indexMatrix =
        Dict{String,
             Vector{UnitRange{Int64}}}(
                 name => UnitRange{Int64}[] for name in keys(jst.children))
    
    for (name, child) in jst.children

        empty!(vector)
        if (!isnothing(child.deltaMap))
            seq = flatten(jst, name)
            vector = BioSequences.findall(query, seq)
        end

        indexMatrix[name] = append!(indexMatrix[name], vector)
        
    end
    return indexMatrix
end

"""
Perform exact search in a Journaled String Tree 2(JST).

# Args: 
   - `jst`: The JSTree2 to search in. 
   - `needle`: The LongDNA{4} sequence to find.

# Returns: 
   - A dictionary mapping node names to exact match ranges.
"""
function exact_search(tree::JSTree2, needle::LongDNA{4})
    # prepare a dict of result‐vectors for each sequence index
    results = Dict(i => UnitRange{Int64}[] for i in 1:tree.length)
    # build the exact‐search query once
    query = ExactSearchQuery(needle)

    # for each “leaf” sequence index...
    for i in 1:tree.length
        # start from the root sequence
        seq = copy(tree.root)

        # apply every delta in ascending position order
        for (pos, buckets) in tree.journal
            for bucket in buckets
                # if this bucket has a mutation for index i, apply it
                if haskey(bucket, Int64(i))
                    entry = bucket[Int64(i)]
                    seq = apply_delta(seq, entry)
                end
            end
        end

        # now search the fully reconstructed seq for exact matches
        matches = BioSequences.findall(query, seq)
        append!(results[i], matches)
    end

    return results
end

function approximate_search(tree::JSTree2, needle::LongDNA{4})
    tol         = ceil(Int, length(needle) * 0.05)
    query       = ApproximateSearchQuery(needle)

    # initial root matches
    root_hits   = approximate_findall(query, tol, tree.root)
    results     = Dict(i => copy(root_hits) for i in 1:tree.length)
    results_set = Dict(i => Set(root_hits)  for i in 1:tree.length)

    for i in 1:tree.length
        applied = false
        seq     = nothing

        for (_pos, buckets) in tree.journal
            for bucket in buckets
                if haskey(bucket, Int64(i))
                    entry      = bucket[Int64(i)]
                    # find which current hits this entry invalidates
                    hit_ranges = collect(results_set[i])
                    invalid    = [r for r in hit_ranges if entry.position in r]

                    # only if there are invalidated root‐hits do we proceed
                    if !isempty(invalid)
                        if !applied
                            seq              = copy(tree.root)
                            results_set[i]   = Set()   # drop all old hits
                            applied          = true
                        end

                        # remove just the invalidated ranges
                        for r in invalid
                            delete!(results_set[i], r)
                        end

                        # apply mutation and add new approximate hits
                        seq = apply_delta(seq, entry)
                        for newr in approximate_findall(query, tol, seq)
                            push!(results_set[i], newr)
                        end
                    end
                end
            end
        end

        # finalize
        results[i] = collect(results_set[i])
    end

    return results
end

function approximate_search(tree::JSTree2, needle::LongDNA{4}, tol::Int64)

    if tol <= 0 || tol >= 100
        error("Tolerance cannot less or 0% or more than 100%")
    end
    tolerance = ceil(Int64, (length(needle) / 100) * tol) 
    query       = ApproximateSearchQuery(needle)

    # initial root matches
    root_hits   = approximate_findall(query, tolerance, tree.root)
    results     = Dict(i => copy(root_hits) for i in 1:tree.length)
    results_set = Dict(i => Set(root_hits)  for i in 1:tree.length)

    for i in 1:tree.length
        applied = false
        seq     = nothing

        for (_pos, buckets) in tree.journal
            for bucket in buckets
                if haskey(bucket, Int64(i))
                    entry      = bucket[Int64(i)]
                    # find which current hits this entry invalidates
                    hit_ranges = collect(results_set[i])
                    invalid    = [r for r in hit_ranges if entry.position in r]

                    # only if there are invalidated root‐hits do we proceed
                    if !isempty(invalid)
                        if !applied
                            seq              = copy(tree.root)
                            results_set[i]   = Set()   # drop all old hits
                            applied          = true
                        end

                        # remove just the invalidated ranges
                        for r in invalid
                            delete!(results_set[i], r)
                        end

                        # apply mutation and add new approximate hits
                        seq = apply_delta(seq, entry)
                        for newr in approximate_findall(query, tol, seq)
                            push!(results_set[i], newr)
                        end
                    end
                end
            end
        end

        # finalize
        results[i] = collect(results_set[i])
    end

    return results
end
### algorithms.jl ends here.
