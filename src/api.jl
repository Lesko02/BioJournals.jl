### -*- Mode: Julia -*-

### api.jl


"""
    sequence(journal, i) -> LongSequence

Materialize derived sequence `i`. This method is part of `FASTX.sequence`,
so a journal and a FASTA/FASTQ record share one function. The result is an
ordinary `LongSequence`: BioSequences `findall`, `translate`, and FASTX
writers accept it directly.

`i` is 1-based. For a `JournaledString` it indexes `deltaMap`. For a `JSTree`
it indexes the sequences tracked by the tree.
"""
function sequence(js :: JournaledString, i :: Integer)
    return apply_delta(js.reference, edits_at(js, i))
end

function sequence(tree :: JSTree, i :: Integer)
    return apply_delta(tree.root, edits_at(tree, i))
end

"""
    sequences(journal) -> Vector{<:LongSequence}

Materialize every derived sequence in `journal`.
"""
function sequences(js :: JournaledString)
    return [sequence(js, i) for i in eachindex(js.deltaMap)]
end

function sequences(tree :: JSTree)
    return [sequence(tree, i) for i in 1:tree.length]
end

"""
    nsequences(journal) -> Int

Number of derived sequences stored in `journal`.
"""
nsequences(js :: JournaledString) = length(js.deltaMap)
nsequences(tree :: JSTree) = tree.length

"""
    get_mutation_history(delta_map) -> String

List every edit in `delta_map` with its timestamp.
"""
function get_mutation_history(delta_map :: DeltaMap)
    io = IOBuffer()
    for (time, entry) in delta_map
        println(io, "Time ", UInt64(time), ": ", entry)
    end
    return String(take!(io))
end


"""
Retrieve mutations within a specific timestamp interval.

# Args: 
   - `delta_map`: The DeltaMap containing mutations. 
   - `time1`: Start timestamp of the interval. 
   - `time2`: End timestamp of the interval.

# Returns: 
   - A string listing mutations within the given range.
"""
function get_mutation_interval(js :: JournaledString,
                               time1 :: Int,
                               time2 :: Int)
    
    io = IOBuffer()
    for i in 1:length(js.deltaMap)
        for (time, entry) in js.deltaMap[i]
            if time1 <= time <= time2
                println(io, "Time ", time, ": ", entry, " at index ", i)
            end
        end
    end
    return String(take!(io))
end


"""
Get sequences at a specific timestamp in a JournaledString.

# Args: 
   - `jst`: The JournaledString to extract sequences from. 
   - `time`: The timestamp limit for reconstruction.

# Returns:
   - A vector of sequences at the given time, in the reference alphabet.
"""
function get_sequences_at_time(jst :: JournaledString{S}, time :: Integer) where {S <: LongSequence}
    limit = Timestamp(time)
    sequences_at_time = Vector{S}(undef, length(jst.deltaMap))
    for i in eachindex(jst.deltaMap)
        kept = JournalEntry{S}[]
        for (entry_time, entry) in jst.deltaMap[i]
            entry_time > limit && break
            push!(kept, entry)
        end
        sequences_at_time[i] = apply_delta(jst.reference, kept)
    end
    return sequences_at_time
end


"""
Simulate a mutation in a JournaledString.

# Args: 
   - `jst`: The JournaledString to modify. 
   - `index`: The index in the delta map. 
   - `entry`: The JournalEntry mutation to add.

# Modifies: 
    The JournaledString in place.
"""
function simulate_mutation!(jst :: JournaledString,
                            index :: Integer,
                            entry :: JournalEntry)
   add_delta!(jst, [index], entry)
end


"""
Remove a mutation at a specific time.

# Args: 
   - `jst`: The JournaledString to modify. 
   - `indices`: vector of indices
   - `time`: The timestamp of the mutation to remove.

# Modifies: 
   - The JournaledString in place.
"""
function remove_mutation!(jst :: JournaledString,
                          index :: Integer,
                          time :: Integer)
    remove_delta!(jst ,index, time)
end


"""
Check if two JournaledStrings are equal.

# Args: 
   - `jst1`: The first JournaledString. 
   - `jst2`: The second JournaledString.

# Returns: 
   - `true` if both have identical reference, time, and deltas.
   - `false` otherwise
"""
function is_equal(jst1 :: JournaledString,
                  jst2 :: JournaledString) :: Bool
    return jst1.reference == jst2.reference &&
           jst1.current_time == jst2.current_time &&
           jst1.reference_id == jst2.reference_id &&
           jst1.deltaMap == jst2.deltaMap
end


"""
Check if two Journaled String Trees (JSTrees) are equal.

`rootChildren` is a cache built by `build_tree!` and is not compared.
Two trees are equal when their reference, clock, length, and journal match.

Args: 
    jst1: The first JSTree. 
    jst2: The second JSTree.

Returns: 
    `true` if both trees describe the same sequences.
"""
function is_equal(jst1 :: JSTree, jst2 :: JSTree) :: Bool
    return jst1.root == jst2.root &&
           jst1.current_time == jst2.current_time &&
           jst1.length == jst2.length &&
           jst1.journal == jst2.journal
end



"""
Print search results from a JournaledString.

# Args: 
   - `results`: A dictionary mapping delta map indices to match ranges.

# Output: 
   - Prints matches or indicates no matches.
"""
function print_results(results :: Dict{Int64, Vector{UnitRange{Int64}}})
    for i in sort!(collect(keys(results)))
        ranges = results[i]
        if isempty(ranges)
            println("No Match in DeltaMap $i")
        else
            println("Match in $i:")
            println("Ranges: ", ranges)
        end
    end
    return nothing
end


"""
Print search results from a Journaled String Tree (JST).

# Args: 
   - `results`: A dictionary mapping node names to match ranges.

# Output: 
   - Prints matches or indicates no matches.
"""
function print_results(results :: Dict{String, Vector{UnitRange{Int64}}})
    for name in sort!(collect(keys(results)))
        ranges = results[name]
        if isempty(ranges)
            println("No Match at $name")
        else
            println("Match at $name:")
            println("Ranges: ", ranges)
        end
    end
    return nothing
end


### api.jl ends here.
