module BioJournals

using BioSequences
using DataStructures
using FASTX

# Definition of DeltaTypes
@enum DeltaType DeltaTypeDel DeltaTypeIns DeltaTypeSnp DeltaTypeSV DeltaTypeCNV

# Custom insert for sequences
function insert!(seq::LongDNA{4}, pos::Int, subseq::LongDNA)
    for symbol in subseq
        Base.insert!(seq, pos, symbol)  # Inserimento dei singoli simboli
    end
    return seq
end

# Custom delete_at for ranges
function delete_at!(seq::LongDNA, pos_range::UnitRange{Int})
    for i in reverse(pos_range)
        Base.deleteat!(seq, i)
    end
    return seq
end

# Naive implementation of structure_variation //TO FIX//
function structure_variation!(seq::LongDNA{4}, pos::Int, subseq::LongDNA)
    # Check if it's a base case of an append!
    if pos == lastindex(seq)
        append!(seq, subseq)
        return seq
    # Check if it's a base case of insert!
    elseif pos < lastindex(seq)
        insert!(seq, pos, subseq)
        return seq
    # Check of an out of bounds insert --> SV
    elseif pos > lastindex(seq)
        i = lastindex(seq)
        range = lastindex(seq):pos-1
        for i in range
            push!(seq, '-')
        end
        append!(seq, subseq)
        return seq
    end
end

# Implementation of CNV
function copy_number_variation!(seq::LongDNA{4}, pos::Int,
    params::Tuple{LongDNA{4}, Int})
   subseq, rep = params  # Deconstruct the tuple
   for i in 1:rep  # Insert subseq `rep` times
       insert!(seq, pos, subseq)
   end    
   return seq
end

# Journal Entry
struct JournalEntry
    delta_type::DeltaType    # Type of delta
    position::Int            # Index of delta
    data::Any                # Additional Data
    time::Int                # Timestamp
end                          

# Definition of JournaledString Structure
mutable struct JournaledString
    reference::LongDNA{4}
    deltaMap::Vector{SortedDict{Int, JournalEntry, Base.Order.ForwardOrdering}}
    current_time::Int
end

# Constructor for JournaledString
function JournaledString(reference::LongDNA{4},
    deltaMap::Vector{SortedDict{Int, JournalEntry, Base.Order.ForwardOrdering}})
    JournaledString(reference, deltaMap, 0)  # Default value for current_time
end

# Function to add a new Delta
function add_delta!(js::JournaledString, indices::Vector{Int}, 
                    delta_type::DeltaType, position::Int, data::Any)
    for idx in indices
       # Create the new JournalEntry
       new_entry = JournalEntry(delta_type, position, data, js.current_time)
        
       # Insert the entry into the SortedDict with `time` as the key
       js.deltaMap[idx][js.current_time] = new_entry

       # Increment the current time
       js.current_time += 1
    end
end

function add_delta!(js::JournaledString,
     indices::Vector{Int}, entry::JournalEntry)
        for idx in indices
        # Insert the entry into the SortedDict with `time` as the key
        js.deltaMap[idx][js.current_time] = entry

        # Increment the current time
        js.current_time += 1
        end
end

function remove_delta!(js::JournaledString, time::Int)
    for idx in deltaMap
        for entry in js.deltaMap[idx]
            if(entry[time] == time)
                delete!(js.deltaMap[idx], time)
            else
            error("No mutation found at time: $time" )
            end
        end
    end
    
end

function apply_delta(reference::LongDNA{4}, 
                    delta::SortedDict{Int, JournalEntry})
    seq = copy(reference)
    for (_, entry) in delta  # Access entries in sorted order of `time`
        # Check on the DeltaType
        if entry.delta_type == DeltaTypeDel
            seq = delete_at!(seq, entry.position:(entry.position + 
                   entry.data - 1))  # Data is the bound of the range
        elseif entry.delta_type == DeltaTypeIns
            seq = insert!(seq, entry.position, LongDNA{4}(entry.data))
             # Single nucleotide permutation
        elseif entry.delta_type == DeltaTypeSnp
            seq[entry.position] = convert(DNA, entry.data)  
            # Larger Structure change
        elseif entry.delta_type == DeltaTypeSV
            seq = structure_variation!(seq, entry.position, entry.data)
        elseif entry.delta_type == DeltaTypeCNV
            seq = copy_number_variation!(seq, entry.position, entry.data)
        end
    end
    return seq
end

# Function to print the n sequences
function print_sequences(jst::JournaledString)
    for i in 1:length(jst.deltaMap)
        modified_seq = apply_delta(jst.reference, jst.deltaMap[i])
        println("Sequence $i: ", modified_seq)
    end
end

# Function to build the n sequences into a sigle string
function build_sequences(jst::JournaledString)
    builded = ""
    for i in 1:length(jst.deltaMap)
        modified_seq = apply_delta(jst.reference, jst.deltaMap[i])
        builded *= "Sequence $i: " * string(modified_seq) * "\n"
    end
    return builded
end

# Function to print all of the deltas
function print_deltas(jst::JournaledString)
    for j in 1:length(jst.deltaMap)
        println("Stringa indice $j:")
        for (time, entry) in jst.deltaMap[j]
            println("  [time=$time] JournalEntry: $entry")
        end
    end
end

#API definitions

function get_mutation_history(delta_map::SortedDict{Int, JournalEntry})
    mutation_history = ""
    for (time, entry) in delta_map
        mutation_history *= "Time $time: $entry\n"
    end
    return mutation_history
end

function get_mutation_interval(delta_map::SortedDict{Int, JournalEntry},
                                 time1::Int, time2::Int)
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
        filtered_delta = SortedDict{Int, JournalEntry}()
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

function is_equal(entry1::JournalEntry, entry2::JournalEntry)::Bool
    return entry1.delta_type == entry2.delta_type &&
           entry1.position == entry2.position &&
           entry1.data == entry2.data
end

function is_equal(jst1::JournaledString, jst2::JournaledString)::Bool

    l1 = length(jst1.deltaMap)
    l2 = length(jst2.deltaMap)

    if js1.reference != js2.reference
       return false
    end
    
    if l1 != l2
        return false
    end

    # Compare each SortedDict in deltaMap
    for (i, dict1) in enumerate(js1.deltaMap)
        dict2 = js2.deltaMap[i]

        # Compare sizes of the SortedDicts
        if length(dict1) != length(dict2)
            return false
        end

        # Compare JournalEntry values, ignoring time
        entries1 = collect(values(dict1))
        entries2 = collect(values(dict2))

        # Ensure the two sets of entries match
        if length(entries1) != length(entries2) ||
           !all(is_equal(entries1[j], entries2[j]) for j in 1:length(entries1))
            return false
        end
    end

    return true

end

export insert!, delete_at!, structure_variation!, copy_number_variation!,
 JournalEntry, JournaledString, add_delta!, remove_delta!, apply_delta,
 print_sequences, build_sequences, print_deltas, get_mutation_history,
 get_mutation_interval, get_sequences_at_time, simulate_mutation!,
 remove_mutation!, is_equal
end
