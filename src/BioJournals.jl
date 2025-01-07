module BioJournals

using BioSequences
using DataStructures
using FASTX

# Definition of DeltaTypes
@enum DeltaType DeltaTypeDel DeltaTypeIns DeltaTypeSnp DeltaTypeSV

# Definition of current_time
current_time = Ref(0)

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
# Journal Entry
struct JournalEntry
    delta_type::DeltaType    # Type of delta
    position::Int            # Index of delta
    data::Any                # Additional Data
    time::Int                # Timestamp
end                          

# Definition of JournaledString Structure
struct JournaledString
    reference::LongDNA{4}                            
    deltaMap::Vector{SortedDict{Int, JournalEntry}}  # Sorted by time
end                                        

# Function to add a new Delta
function add_delta!(deltaMap, indices::Vector{Int}, 
                    delta_type::DeltaType, position::Int, data::Any)
    for idx in indices
       # Create the new JournalEntry
       new_entry = JournalEntry(delta_type, position, data, current_time[])
        
       # Insert the entry into the SortedDict with `time` as the key
       deltaMap[idx][current_time[]] = new_entry

       # Increment the current time
       current_time[] += 1
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

# Function to print all of the deltas
function print_deltas(jst::JournaledString)
    for j in 1:length(jst.deltaMap)
        println("Stringa indice $j:")
        for (time, entry) in jst.deltaMap[j]
            println("  [time=$time] JournalEntry: $entry")
        end
    end
end


export JournalEntry, JournaledString, add_delta!, apply_delta, print_sequences,
        print_deltas
end
