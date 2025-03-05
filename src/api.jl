"""
Retrieve the full mutation history from a DeltaMap.

Args: 
    delta_map: The DeltaMap containing mutations.

Returns: 
    A string listing all mutations with timestamps.
"""
function get_mutation_history(delta_map::DeltaMap)
    mutation_history = ""
    for (time, entry) in delta_map
        mutation_history *= "Time $time: $entry\n"
    end
    return mutation_history
end

"""
Retrieve mutations within a specific time interval.

Args: 
    delta_map: The DeltaMap containing mutations. 
    time1: Start time of the interval. 
    time2: End time of the interval.

Returns: 
    A string listing mutations within the given range.
"""
function get_mutation_interval(delta_map::DeltaMap, time1::Int, time2::Int)
    mutation_interval = ""
    for (time, entry) in delta_map
        if time1 <= time <= time2
            mutation_interval *= "Time $time: $entry\n"
        end
    end
    return mutation_interval
end

"""
Get sequences at a specific time in a JournaledString.

Args: 
    jst: The JournaledString to extract sequences from. 
    time: The time point for reconstruction.

Returns: 
    A vector of LongDNA{4} sequences at the given time.
"""
function get_sequences_at_time(jst::JournaledString, time::Int)
    sequences_at_time = Vector{LongDNA{4}}(undef, length(jst.deltaMap))
    for i in 1:length(jst.deltaMap)
        filtered_delta = DeltaMap()
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

Args: 
    jst: The JournaledString to modify. 
    index: The index in the delta map. 
    entry: The JournalEntry mutation to add.

Modifies: 
    The JournaledString in place.
"""
function simulate_mutation!(jst::JournaledString, index::Int,
     entry::JournalEntry)
   add_delta!(jst.deltaMap, index, entry)
end

"""
Remove a mutation at a specific time.

Args: 
    jst: The JournaledString to modify. 
    indices: vector of indices
    time: The timestamp of the mutation to remove.

Modifies: 
    The JournaledString in place.
"""
function remove_mutation!(jst::JournaledString, indices::Vector{Int}, time::Int)
    remove_delta!(jst.deltaMap,indices, time)
end

"""
Check if two JournaledStrings are equal.

Args: 
    jst1: The first JournaledString. 
    jst2: The second JournaledString.

Returns: 
    `true` if both have identical reference, time, and deltas.
"""
function is_equal(jst1::JournaledString, jst2::JournaledString)::Bool
    if hash(jst1.reference) != hash(jst2.reference)
        return false
    end
    if jst1.current_time != jst2.current_time
        return false
    end
return hash(jst1.deltaMap)==hash(jst2.deltaMap)
end

"""
Print search results from a JournaledString.

Args: 
    results: A dictionary mapping delta map indices to match ranges.

Output: 
    Prints matches or indicates no matches.
"""
function print_results(results::Dict{Int64, Vector{UnitRange{Int64}}})
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

Args: 
    results: A dictionary mapping node names to match ranges.

Output: 
    Prints matches or indicates no matches.
"""
function print_results(results::Dict{String, Vector{UnitRange{Int64}}})
    for (name, _) in results
        if isempty(results[name])
            println("No Match at $name")
        else
            println("Match at $name:")
            println("Ranges: ", results[name])
        end
    end
end