### -*- Mode: Julia -*-

### functions.jl


"""
Insert a subsequence into a DNA sequence.
Uses the BioJournals insert.

# Args:
   - `seq`: A LongDNA{4} sequence to modify. 
   - `pos`: Position where the subsequence is inserted. 
   - `subseq`: A LongDNA sequence to insert.

# Returns:
   - Updated LongDNA{4} sequence.

# Examples:
```julia-repl
    julia> insert!(LongDNA{4}("ACGT"), 2, LongDNA("CG")) 
    LongDNA{4}("ACCGT") 
```
"""
function insert!(seq :: LongDNA{4}, pos :: Int, subseq :: LongDNA)
    for symbol in reverse(subseq)
        BioSequences.insert!(seq, pos, symbol)
    end
    return seq
end


"""
Delete a range of symbols from a DNA sequence.

# Args:
   - `seq`: A LongDNA sequence to modify. 
   - `pos_range`: Range of positions to delete.

# Returns:
   - Updated LongDNA sequence.

# Examples:
```julia-repl
    julia> delete_at!(LongDNA("ACGTACGT"), 3:5) 
    LongDNA("ACACGT") 
```
"""
function delete_at!(seq :: LongDNA, pos_range :: UnitRange{Int})
    for i in reverse(pos_range)
        BioSequences.deleteat!(seq, i)
    end
    return seq
end


"""
Apply a structure variation to a DNA sequence.

# Args: 
   - `seq`: A LongDNA{4} sequence to modify. 
   - `pos`: Position where the variation occurs. 
   - `subseq`: LongDNA sequence for variation.

# Returns: 
   - Updated LongDNA{4} sequence.

# Examples: 
```julia-repl
    julia> structure_variation!(LongDNA{4}("ACGT"), 4, LongDNA("TT")) 
    LongDNA{4}("ACGTTT") 
```
"""
function structure_variation!(seq :: LongDNA{4},
                              pos :: Int,
                              subseq :: LongDNA)
    ## Check if it's a base case of an append!
    if pos == lastindex(seq)
        append!(seq, subseq)
        return seq
        ## Check if it's a base case of insert!
    elseif pos < lastindex(seq)
        insert!(seq, pos, subseq)
        return seq
        ## Check of an out of bounds insert --> SV
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

# Args: 
   - `seq`: A LongDNA{4} sequence to modify. 
   - `pos`: Position for the variation. 
   - `params`: Tuple (subseq, repetitions).

# Returns: 
   - Updated LongDNA{4} sequence.

# Examples: 
```julia-repl
    julia> copy_number_variation!(LongDNA{4}("ACGT"), 2, (LongDNA{4}("GG"), 3)) 
    LongDNA{4}("ACGGGGGGT") 
```    
"""
function copy_number_variation!(seq :: LongDNA{4},
                                pos :: Int,
                                params :: Tuple{LongDNA{4}, Int})
    subseq, rep = params  # Deconstruct the tuple
    for i in 1:rep  # Insert subseq `rep` times
        insert!(seq, pos, subseq)
    end    
    return seq
end


### End of deltas


"""
Flattens a JSTree to obtain a sequence.

# Args: 
   - `tree`: JSTree structure. 
   - `node_name`: Node to flatten.

# Returns: 
   - Flattened sequence.

# Examples: 
```julia-repl
    julia> flatten(tree, "child_node") 
    LongDNA{4}("ACGT") 
```
"""
function flatten(tree :: JSTree, node_name :: String)
    ## Base case: Root sequence
    if node_name == "root"
        return tree.root
    end

    ## Recursive case: Flatten parent and apply deltas
    node = tree.children[node_name]
    parent_sequence = flatten(tree, node.parent.name)
    return apply_delta(parent_sequence, node.deltaMap)
end


"""
Prints a visual tree representation.

# Args: 
   - `tree`: JSTree structure. 
   - `node_name`: Node to start from. 
   - `indent`: Indentation level.

# Returns: 
   - None (prints tree).

# Examples: 
```julai-repl
    julia> print_tree(tree) 
    |- root
       |- child1
       |- child2 
```
"""
function print_tree(tree :: JSTree,
                    node_name :: String = "root",
                    indent :: Int = 0)
    println(" "^(indent * 2) * "|- " * node_name) 

    for (child_name, child_node) in tree.children
        
        if child_node.parent !== nothing && child_node.parent.name == node_name
            print_tree(tree, child_name, indent + 1) 
        end
    end
end

function print_tree2(tree::JSTree2)

    println("|- root")
    
    for (pos, node) in sort(collect(pairs(tree.children)), by = x -> x.first)
        print_position_node(pos, node, 1)
    end
end

function print_position_node(pos::Int, node::JSTNode2, indent::Int)

    indent_str = "  " * "  "^(indent - 1)
    

    println(indent_str * "|- pos_$pos", collect(keys(node.deltaMap)), 
                                        collect(values(node.deltaMap))) 
    # Recursively print children of the current node
    for (child_pos, child_node) in sort(collect(pairs(node.children)),
         by = x -> x.first)
        print_position_node(child_pos, child_node, indent + 1)
    end
end

function reconstruct(tree::JSTree2)
    results = Dict{Int, LongDNA{4}}()
    results[0] = tree.root
    counter = 1;

    function explore(node::JSTNode2, current_deltas::DeltaMap)

        for d in node.deltaMap
            push!(current_deltas, d)
        end

        if !isempty(node.children)
            results[counter] = apply_delta(tree.root, current_deltas)
            counter += 1
            for (_, child) in node.children
                explore(child, current_deltas)
            end
        end
    end

    # Parte dai nodi direttamente collegati alla root
    for (_, node) in tree.children
        explore(node, DeltaMap())
    end

    return results

end
"""
Prints all JSTree sequences.

# Args:
   - `tree`: JSTree structure.

# Returns:
   - None (prints sequences).

# Examples:
```julia-repl
    julia> print_sequences(tree) 
    root: ACGT... 
    child1: ACG... 
```
"""
function print_sequences(tree :: JSTree)
    println("root: ", tree.root)
    for (name, node) in tree.children
        if name != "root"
            println("$name: ")
            println(flatten(tree, node.name))
        end
    end
end


"""
Prints sequences from a JournaledString.

# Args: 
   - `jst`: JournaledString structure.

# Returns: 
   - None (prints sequences).

# Examples: 
```julia-repl
    julia> print_sequences(jst) 
    Sequence 1: ACGT... 
```
"""
function print_sequences(jst::JournaledString)
    for i in 1:length(jst.deltaMap)
        modified_seq = apply_delta(jst.reference, jst.deltaMap[i])
        println("Sequence $i: ", modified_seq)
    end
end


"""
Builds a string with all sequences.

# Args:
   - `js`: JournaledString structure.

# Returns:
   - Concatenated sequences.

# Examples:
```julia-repl
    julia> build_sequences(js) 
    "Sequence 1: ACGT...
    Sequence 2: TGCA...
    ..."
```    
"""
function build_sequences(js :: JournaledString)
    builded = ""
    for i in 1:length(js.deltaMap)
        modified_seq = apply_delta(js.reference, js.deltaMap[i])
        builded *= "Sequence $i: " * string(modified_seq) * "\n"
    end
    return builded
end


"""
Prints all deltas in a JournaledString.

# Args:
   - `js`: JournaledString structure.

# Returns:
   - None (prints deltas).

# Examples:
```julia-repl
    julia> print_deltas(js)
    " String number 1:
      [time=1] JournalEntry: ... "
```
"""
function print_deltas(js :: JournaledString)
    for j in 1:length(js.deltaMap)
        println("String number $j:")
        for (time, entry) in js.deltaMap[j]
            time_int = UInt64(time)  # Explicit conversion
            println("[time = $time_int] JournalEntry: $entry")
        end
    end
end


"""
Adds a delta to a JournaledString.

# Args:
   - `js`: JournaledString structure. 
   - `indices`: Indices to update. 
   - `delta_type`: Type of delta. 
   - `position`: Position for delta. 
   - `data`: Data for the delta.

# Returns:
   - None (modifies JournaledString).

# Examples:
```julia-repl
    julia> add_delta!(js, [1, 2], DeltaTypeIns, 3, "AC") 
```    
"""
function add_delta!(js :: JournaledString,
                    indices :: Vector{Int}, 
                    delta_type :: DeltaType,
                    position :: Int,
                    data :: Any)
    for idx in indices
        ## Create the new JournalEntry
        new_entry = JournalEntry(delta_type, position, data, js.current_time)

        ## Insert the entry into the SortedDict with `time` as the key
        js.deltaMap[idx][js.current_time] = new_entry

        ## Increment the current time
        js.current_time += UInt64(1)
    end
end

#new experimental for jst2
function add_delta!(tree :: JSTree2,
                    indices :: Vector{Int}, 
                    delta_type :: DeltaType,
                    position :: Int,
                    data :: Any)

    for idx in indices
        ## Create the new JournalEntry
        new_entry = JournalEntry(delta_type, position, data, tree.current_time)


        entries = get!(tree.journal, Int64(position)) do
        fill(nothing, tree.length)
        end

        entries[idx] = new_entry
        ## Increment the current time
        tree.current_time += UInt64(1)

    end

    
end


"""
Adds a JournalEntry to indices.

# Args:
   - `js`: JournaledString structure.
   - `indices`: Indices to update.
   - `entry`: JournalEntry object.

# Returns:
   - None (modifies JournaledString).

# Examples:
```julia-repl
    julia> add_delta!(js, [1], entry)
```
"""
function add_delta!(js :: JournaledString,
                    indices :: Vector{Int}, 
                    entry :: JournalEntry)
    for idx in indices
        js.deltaMap[idx][js.current_time] = entry
        js.current_time += UInt64(1)
    end
end


"""
Add a mutation entry (delta) to a node in a Journaled String Tree (JSTree).

# Args:
   - `jst`: The JSTree to modify.
   - `node_name`: The name of the node to which the delta is added.
   - `delta_type`: The type of mutation (e.g., insertion, deletion).
   - `position`: The position in the sequence where the mutation occurs.
   - `data`: The data associated with the mutation.

# Errors:
   - Throws an error if attempting to add a delta to the root node.
   - Throws an error if the specified node does not exist in the JSTree.
   - Throws an error if the specified node lacks a delta map.
"""
function add_delta!(jst :: JSTree,
                    node_name :: String,
                    delta_type::DeltaType,
                    position :: Int,
                    data :: Any)

    if node_name == "root"
        error("Cannot add delta to the root node.")
    end
    
    if !haskey(jst.children, node_name)
        error("Node '$node_name' does not exist.")
    end
    
    node = jst.children[node_name]
    deltaMap = node.deltaMap
    next_time = isempty(deltaMap) ? Timestamp(0) : maximum(keys(deltaMap)) + UInt64(1)
    new_entry = JournalEntry(delta_type, position, data, next_time)
    deltaMap[next_time] = new_entry
end


"""
Add a mutation entry (delta) to a node in a Journaled String Tree (JSTree).

# Args:
   - `jst`: The JSTree to modify.
   - `node_name`: The name of the node to which the delta is added.
   - `entry`: The JournalEntry to add.

# Errors:
   - Throws an error if attempting to add a delta to the root node.
   - Throws an error if the specified node does not exist in the JSTree.
   - Throws an error if the specified node lacks a delta map.
"""
function add_delta!(jst :: JSTree, node_name :: String, entry :: JournalEntry)

    if node_name == "root"
        error("Cannot add delta to the root node.")
    end
    
    if !haskey(jst.children, node_name)
        error("Node '$node_name' does not exist.")
    end
    
    node = jst.children[node_name]
    deltaMap = node.deltaMap
    next_time = isempty(deltaMap) ? Timestamp(0) : maximum(keys(deltaMap)) + UInt64(1)
    new_entry = JournalEntry(entry.delta_type,
                             entry.position,
                             entry.data,
                             next_time)  
    deltaMap[next_time] = new_entry
end


"""
Removes a delta by time key.

# Args: 
   - `js`: JournaledString structure. 
   - `time`: Time key of the delta.

# Returns: 
   - None (modifies JournaledString).

# Examples: 
```julia-repl
    julia> remove_delta!(js, 5) 
```
"""
function remove_delta!(js :: JournaledString, idx :: Int, time :: Int)
    deletioncc = 0
    for (timer , entry) in js.deltaMap[idx]
        if timer == time
            delete!(js.deltaMap[idx], UInt(time))
            deletioncc += 1
        end
    end

    if deletioncc == 0
        error("No mutation found at time: $time at index: $idx" )
    end
end

# needs overloading
function add_delta2!(tree::JSTree2, entry::JournalEntry)
    pos = entry.position
    current_node = tree.children
    path = []

    while pos in keys(current_node)
        path_node = current_node[pos]
        path = [path..., path_node]
        current_node = path_node.children
    end

    # new node at position
    parent = isempty(path) ? nothing : last(path)
    deltaMap = DeltaMap()
    deltaMap[entry.time] = entry
    children = Dict{Int64, JSTNode2}()
    new_node = JSTNode2(parent, deltaMap, children)
    
    current_node[pos] = new_node
end
"""
Removes a delta by time key.

# Args:
   - `jst`: JSTree structure.
   - `node_name`: Name of the node.
   - `time`: Time key of the delta.

# Returns:
   - None (modifies JSTNode).

# Examples:
```julia-repl
    julia> remove_delta!(jst, "children1", 5)
```
"""
function remove_delta!(jst :: JSTree, node_name :: String, time :: Int)
    deletioncc = 0
    if node_name == "root"
        error("Cannot delete delta of the root node.")
    end
    if !haskey(jst.children, node_name)
        error("Node '$node_name' does not exist.")
    end

    for (timer , entry) in jst.children[node_name].deltaMap
        if timer == time
            delete!(jst.children[node_name].deltaMap, UInt(time))
            deletioncc += 1
        end
    end
    if deletioncc == 0
        error("No mutation found at time: $time at child: $node_name" )
    end
end


"""
Apply a delta to a sequence.

# Args:
   - `reference`: LongDNA{4} reference sequence.
   - `entry`: JournalEntry delta.

# Returns:
   - Modified LongDNA{4} sequence.

# Examples:
```julia-repl
    julia> apply_delta(ref_seq, entry)
```
"""
function apply_delta(reference :: LongDNA{4}, entry :: JournalEntry)
    seq = copy(reference)

    ## Check on the DeltaType
    if entry.delta_type == DeltaTypeDel
        seq = delete_at!(seq,
                         entry.position:(entry.position + 
                             entry.data - 1))  # Data is the bound of the range
    elseif entry.delta_type == DeltaTypeIns
        seq = insert!(seq, entry.position, LongDNA{4}(entry.data))

    elseif entry.delta_type == DeltaTypeSnp
        ## Single nucleotide permutationa
        seq[entry.position] = convert(DNA, entry.data)  

    elseif entry.delta_type == DeltaTypeSV
        ## Larger Structure change
        seq = structure_variation!(seq, entry.position, entry.data)
    elseif entry.delta_type == DeltaTypeCNV
        seq = copy_number_variation!(seq, entry.position, entry.data)
    end
    return seq
end


"""
Apply multiple deltas to a sequence.

# Args: 
   - `reference`: LongDNA{4} reference sequence. 
   - `delta`: DeltaMap of JournalEntries.

# Returns: 
   - Modified LongDNA{4} sequence.

# Examples: 
```julia-repl
    julia> apply_delta(ref_seq, delta_map) 
```
"""
function apply_delta(reference :: LongDNA{4}, delta :: DeltaMap)

    seq = copy(reference)
    for (_, entry) in delta 
        ## Check on the DeltaType
        if entry.delta_type == DeltaTypeDel
            seq = delete_at!(seq, entry.position:(entry.position + 
            entry.data - 1))  # Data is the bound of the range
        elseif entry.delta_type == DeltaTypeIns
            seq = insert!(seq, entry.position, LongDNA{4}(entry.data))
        elseif entry.delta_type == DeltaTypeSnp
            ## Single nucleotide permutation

            seq[entry.position] = convert(DNA, entry.data)  

        elseif entry.delta_type == DeltaTypeSV
            ## Larger Structure change
            seq = structure_variation!(seq, entry.position, entry.data)
        elseif entry.delta_type == DeltaTypeCNV
            seq = copy_number_variation!(seq, entry.position, entry.data)
        end
    end
    return seq
end


### functions.jl ends here.
