### -*- Mode: Julia -*-

### api.jl


"""
Retrieve the full mutation history from a DeltaMap.

# Args: 
   - `delta_map`: The DeltaMap containing mutations.

# Returns: 
   - A string listing all mutations with timestamps.
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
   - A vector of LongDNA{4} sequences at the given time.
"""
function get_sequences_at_time(jst :: JournaledString, time :: Int)
    sequences_at_time = Vector{LongDNA{4}}(undef, length(jst.deltaMap))
    for i in 1:length(jst.deltaMap)
        filtered_delta = DeltaMap(Base.Order.Forward)
        for (entry_time, entry) in jst.deltaMap[i]
            if entry_time > time
                break
            end
            filtered_delta[entry_time] = entry
        end
        sequences_at_time[i] = apply_delta(jst.reference, filtered_delta)
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
                            index :: Int,
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
                          index :: Int,
                          time  :: Int)
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
    for i in 1:length(results)
        if isempty(results[i])
            println("No Match in DeltaMap $i")
        else
            println("Match in $i:")
            println("Ranges: ", results[i])
        end
    end
end


"""
Print search results from a Journaled String Tree (JST).

# Args: 
   - `results`: A dictionary mapping node names to match ranges.

# Output: 
   - Prints matches or indicates no matches.
"""
function print_results(results :: Dict{String, Vector{UnitRange{Int64}}})
    for (name, _) in results
        if isempty(results[name])
            println("No Match at $name")
        else
            println("Match at $name:")
            println("Ranges: ", results[name])
        end
    end
end


### api.jl ends here.
