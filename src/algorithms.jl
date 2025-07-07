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
function exact_search(tree::JSTree, needle::LongDNA{4})
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

"""
Perform approximate search in a Journaled String Tree (JST).

# Args:
   - `jst`: The JSTree to search in.
   - `needle`: The LongDNA{4} sequence to find.

# Returns:
    A dictionary mapping node names to match ranges.
"""

function approximate_search(tree::JSTree, needle::LongDNA{4})

    tol = ceil(Int64, (length(needle) / 100) * 5) 
    query = ApproximateSearchQuery(needle)

    # Initial root matches
    root_hits = approximate_findall(query, tol, tree.root)
    results = Dict(i => copy(root_hits) for i in 1:tree.length)
    results_set = Dict(i => Set(root_hits) for i in 1:tree.length)

    # Collect journal positions once (outside sequence loop)
    journal_positions = collect(keys(tree.journal))
    sort!(journal_positions)  # Ensure processing in position order

    for i in 1:tree.length
        applied = false
        current_seq = nothing
        current_positions = copy(journal_positions)  # Positions to process for this sequence

        # Process journal entries in position order
        for pos in current_positions
            buckets = tree.journal[pos]
            for bucket in buckets
                if haskey(bucket, Int64(i))
                    entry = bucket[Int64(i)]
                    # Check if this entry affects existing hits
                    hit_ranges = collect(results_set[i])
                    invalid = [r for r in hit_ranges if entry.position in r]

                    if !isempty(invalid)
                        if !applied
                            # Initialize sequence with root only once
                            current_seq = copy(tree.root)
                            results_set[i] = Set()
                            applied = true
                        end

                        # Remove invalidated ranges
                        for r in invalid
                            delete!(results_set[i], r)
                        end

                        # Apply mutation and get new hits
                        current_seq = apply_delta(current_seq, entry)
                        new_hits = approximate_findall(query, tol, current_seq)
                        
                        # Add new hits and break out after applying mutation
                        for newr in new_hits
                            push!(results_set[i], newr)
                        end
                        
                        # BREAK AFTER APPLYING MUTATION
                        break  # Exit bucket loop after applying mutation
                    end
                end
            end
        end

        # Finalize results for this sequence
        results[i] = sort!(collect(results_set[i]), by=first)
    end

    return results
end

"""
Perform approximate search in a Journaled String Tree (JST) with tolerance.

# Args:
   - `jss`: The JournaledString to search in.
   - `needle`: The LongDNA{4} sequence to find.
   - `tol`: Allowed mismatch percentage (1-99%).

# Returns:
   - A dictionary mapping time points to match ranges.

# Errors:
   - Throws an error if tolerance is ≤0% or ≥100%.
"""
function approximate_search(tree::JSTree, needle::LongDNA{4}, tol::Int64)
    if tol <= 0 || tol >= 100
        error("Tolerance cannot less or 0% or more than 100%")
    end
    tolerance = ceil(Int64, (length(needle) / 100) * tol) 
    query = ApproximateSearchQuery(needle)

    # Initial root matches
    root_hits = approximate_findall(query, tolerance, tree.root)
    results = Dict(i => copy(root_hits) for i in 1:tree.length)
    results_set = Dict(i => Set(root_hits) for i in 1:tree.length)

    # Collect journal positions once (outside sequence loop)
    journal_positions = collect(keys(tree.journal))
    sort!(journal_positions)  # Ensure processing in position order

    for i in 1:tree.length
        applied = false
        current_seq = nothing
        current_positions = copy(journal_positions)  # Positions to process for this sequence

        # Process journal entries in position order
        for pos in current_positions
            buckets = tree.journal[pos]
            for bucket in buckets
                if haskey(bucket, Int64(i))
                    entry = bucket[Int64(i)]
                    # Check if this entry affects existing hits
                    hit_ranges = collect(results_set[i])
                    invalid = [r for r in hit_ranges if entry.position in r]

                    if !isempty(invalid)
                        if !applied
                            # Initialize sequence with root only once
                            current_seq = copy(tree.root)
                            results_set[i] = Set()
                            applied = true
                        end

                        # Remove invalidated ranges
                        for r in invalid
                            delete!(results_set[i], r)
                        end

                        # Apply mutation and get new hits
                        current_seq = apply_delta(current_seq, entry)
                        new_hits = approximate_findall(query, tolerance, current_seq)
                        
                        # Add new hits and break out after applying mutation
                        for newr in new_hits
                            push!(results_set[i], newr)
                        end
                        
                        # BREAK AFTER APPLYING MUTATION
                        break  # Exit bucket loop after applying mutation
                    end
                end
            end
        end

        # Finalize results for this sequence
        results[i] = sort!(collect(results_set[i]), by=first)
    end

    return results
end

### algorithms.jl ends here.
