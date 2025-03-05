"""
Insert a subsequence into a DNA sequence.
Uses the BioJournals insert.

Args: 
    seq: A LongDNA{4} sequence to modify. 
    pos: Position where the subsequence is inserted. 
    subseq: A LongDNA sequence to insert.

Returns: Updated LongDNA{4} sequence.

Examples: 
    insert!(LongDNA{4}("ACGT"), 2, LongDNA("CG")) LongDNA{4}("ACCGT") 
"""

function insert!(seq::LongDNA{4}, pos::Int, subseq::LongDNA)
    for symbol in subseq
        BioJournals.insert!(seq, pos, symbol)
    end
    return seq
end

"""
Delete a range of symbols from a DNA sequence.

Args: 
    seq: A LongDNA sequence to modify. 
    pos_range: Range of positions to delete.

Returns: 
    Updated LongDNA sequence.

Examples: 
    delete_at!(LongDNA("ACGTACGT"), 3:5) 
    LongDNA("ACACGT") 
"""
function delete_at!(seq::LongDNA, pos_range::UnitRange{Int})
    for i in reverse(pos_range)
        BioJournals.deleteat!(seq, i)
    end
    return seq
end

"""
Apply a structure variation to a DNA sequence.

Args: 
    seq: A LongDNA{4} sequence to modify. 
    pos: Position where the variation occurs. 
    subseq: LongDNA sequence for variation.

Returns: 
    Updated LongDNA{4} sequence.

Examples: 
    structure_variation!(LongDNA{4}("ACGT"), 4, LongDNA("TT")) 
    LongDNA{4}("ACGTTT") 
"""

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

"""
Apply a copy number variation (CNV).

Args: 
    seq: A LongDNA{4} sequence to modify. 
    pos: Position for the variation. 
    params: Tuple (subseq, repetitions).

Returns: 
    Updated LongDNA{4} sequence.

Examples: 
    copy_number_variation!(LongDNA{4}("ACGT"), 2, (LongDNA{4}("GG"), 3)) 
    LongDNA{4}("ACGGGGGGT") 
"""
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
"""
Flatten a JSTree to obtain a sequence.

Args: 
    tree: JSTree structure. 
    node_name: Node to flatten.

Returns: 
    Flattened sequence.

Examples: 
    flatten(tree, "child_node") 
    LongDNA{4}("ACGT") 
"""

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

"""
Print a visual tree representation.

Args: 
    tree: JSTree structure. 
    node_name: Node to start from. 
    indent: Indentation level.

Returns: 
    None (prints tree).

Examples: 
    print_tree(tree) 
    |- root
       |- child1
       |- child2 
"""
function print_tree(tree::JSTree, node_name::String = "root", indent::Int = 0)
    println(" "^(indent * 2) * "|- " * node_name) 

    for (child_name, child_node) in tree.children
        
        if child_node.parent !== nothing && child_node.parent.name == node_name
            print_tree(tree, child_name, indent + 1) 
        end
    end
end

"""
Print all JSTree sequences.

Args: 
    tree: JSTree structure.

Returns: 
    None (prints sequences).

Examples: 
    print_sequences(tree) 
    root: ACGT... 
    child1: ACG... 
"""
function print_sequences(tree::JSTree)
    println("root: ", tree.root)
    for (name, node) in tree.children
        if name != "root"
        println("$name: ")
        println(flatten(tree, node.name))
        end
    end
end

"""
Print sequences from a JournaledString.

Args: 
    jst: JournaledString structure.

Returns: 
    None (prints sequences).

Examples: 
    print_sequences(jst) 
    Sequence 1: ACGT... 
"""

function print_sequences(jst::JournaledString)
    for i in 1:length(jst.deltaMap)
        modified_seq = apply_delta(jst.reference, jst.deltaMap[i])
        println("Sequence $i: ", modified_seq)
    end
end

"""
Build a string with all sequences.

Args: 
    js: JournaledString structure.

Returns: 
    Concatenated sequences.

Examples: 
    build_sequences(js) 
    "Sequence 1: ACGT...
    Sequence 2: TGCA...
    " 
"""
function build_sequences(js::JournaledString)
    builded = ""
    for i in 1:length(js.deltaMap)
        modified_seq = apply_delta(js.reference, js.deltaMap[i])
        builded *= "Sequence $i: " * string(modified_seq) * "\n"
    end
    return builded
end

"""
Print all deltas in a JournaledString.

Args: 
    js: JournaledString structure.

Returns: 
    None (prints deltas).

Examples: 
    print_deltas(js) 
    String number 1: 
      [time=1] JournalEntry: ... 
"""
function print_deltas(js::JournaledString)
    for j in 1:length(js.deltaMap)
        println("String number $j:")
        for (time, entry) in js.deltaMap[j]
            println("  [time=$time] JournalEntry: $entry")
        end
    end
end

"""
Add a delta to a JournaledString.

Args: 
    js: JournaledString structure. 
    indices: Indices to update. 
    delta_type: Type of delta. 
    position: Position for delta. 
    data: Data for the delta.

Returns: 
    None (modifies JournaledString).

Examples: 
    add_delta!(js, [1, 2], DeltaTypeIns, 3, "AC") 
"""
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

"""
Add a JournalEntry to indices.

Args: 
    js: JournaledString structure. 
    indices: Indices to update. 
    entry: JournalEntry object.

Returns: 
    None (modifies JournaledString).

Examples: 
    add_delta!(js, [1], entry) 
"""
function add_delta!(js::JournaledString, indices::Vector{Int}, 
    entry::JournalEntry)
for idx in indices

    js.deltaMap[idx][js.current_time] = entry

    js.current_time += 1
    end
end

"""
Remove a delta by time key.

Args: 
    js: JournaledString structure. 
    time: Time key of the delta.

Returns: 
    None (modifies JournaledString).

Examples: 
    remove_delta!(js, 5) 
"""
function remove_delta!(js::JournaledString, indices::Vector{Int}, time::Int)
    for idx in indices
        for entry in js.deltaMap[idx]
            if entry[time] == time
                delete!(js.deltaMap[idx], time)
            else
                error("No mutation found at time: $time" )
            end
        end
    end
end

"""
Apply a delta to a sequence.

Args: 
    reference: LongDNA{4} reference sequence. 
    entry: JournalEntry delta.

Returns: 
    Modified LongDNA{4} sequence.

Examples: 
    apply_delta(ref_seq, entry) 
"""
function apply_delta(reference::LongDNA{4}, entry::JournalEntry)
    seq = copy(reference)
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
    return seq
end

"""
Apply multiple deltas to a sequence.

Args: 
    reference: LongDNA{4} reference sequence. 
    delta: DeltaMap of JournalEntries.

Returns: 
    Modified LongDNA{4} sequence.

Examples: 
    apply_delta(ref_seq, delta_map) 
"""
function apply_delta(reference::LongDNA{4}, delta::DeltaMap)

    seq = copy(reference)
    for (_, entry) in delta 
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
