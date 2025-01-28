#API definitions

function get_mutation_history(delta_map::DeltaMap)
    mutation_history = ""
    for (time, entry) in delta_map
        mutation_history *= "Time $time: $entry\n"
    end
    return mutation_history
end

function get_mutation_interval(delta_map::DeltaMap, time1::Int, time2::Int)
    mutation_interval = ""
    for (time, entry) in delta_map
        if time1 <= time <= time2
            mutation_interval *= "Time $time: $entry\n"
        end
    end
    return mutation_interval
end

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

function simulate_mutation!(jst::JournaledString, index::Int,
     entry::JournalEntry)
   add_delta!(jst.deltaMap, index, entry)
end

function remove_mutation!(jst::JournaledString, time::Int)
    remove_delta!(jst.deltaMap, time)
end

function is_equal(jst1::JournaledString, jst2::JournaledString)::Bool
    if hash(jst1.reference) != hash(jst2.reference)
        return false
    end
    if jst1.current_time != jst2.current_time
        return false
    end
return hash(jst1.deltaMap)==hash(jst2.deltaMap)
end