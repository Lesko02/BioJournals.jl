# Custom insert for sequences
function insert!(seq::LongDNA{4}, pos::Int, subseq::LongDNA)
    for symbol in subseq
        Base.insert!(seq, pos, symbol)  # Inserimento dei singoli simboli
    end
    return seq
end

"""
delete_at!

Custom delete_at! for ranges
"""
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

#end of deltas

function flatten(tree::JSTree, node_name::String)
    # Base case: Root sequence
    if node_name == "root"
        return tree.root
    end

    # Recursive case: Flatten parent and apply deltas
    node = tree.children[node_name]
    parent_sequence = flatten(tree, node.parent.name)
    return apply_delta(parent_sequence, node.deltaMap)
end

function print_tree(tree::JSTree, node_name::String = "root", indent::Int = 0)
    println(" "^(indent * 2) * "|- " * node_name)  # Correct string concatenation
    
    # Recursively print children
    for (child_name, child_node) in tree.children
        # Only print children whose parent matches the current node's name
        if child_node.parent !== nothing && child_node.parent.name == node_name
            print_tree(tree, child_name, indent + 1)  # Recursively print children
        end
    end
end

function print_sequences(tree::JSTree)
    println("root: ", tree.root)
    for (name, node) in tree.children
        if name != "root"
        println("$name: ")
        println(flatten(tree, node.name))
        end
    end
end


# Function to print the n sequences
function print_sequences(jst::JournaledString)
    for i in 1:length(jst.deltaMap)
        modified_seq = apply_delta(jst.reference, jst.deltaMap[i])
        println("Sequence $i: ", modified_seq)
    end
end

# Function to build the n sequences into a sigle string
function build_sequences(js::JournaledString)
    builded = ""
    for i in 1:length(js.deltaMap)
        modified_seq = apply_delta(js.reference, js.deltaMap[i])
        builded *= "Sequence $i: " * string(modified_seq) * "\n"
    end
    return builded
end

# Function to print all of the deltas
function print_deltas(js::JournaledString)
    for j in 1:length(js.deltaMap)
        println("Stringa indice $j:")
        for (time, entry) in js.deltaMap[j]
            println("  [time=$time] JournalEntry: $entry")
        end
    end
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
