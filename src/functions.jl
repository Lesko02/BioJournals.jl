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
    julia> String(insert!(LongDNA{4}("ACGT"), 2, LongDNA{4}("CG")))
    "ACGCGT"
```
"""
function insert!(seq :: LongDNA{4}, pos :: Int, subseq :: LongDNA)
    # One shift of the tail, instead of one insert per symbol.
    return BioSequences.spliceinto!(seq, pos, subseq)
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
    julia> String(delete_at!(LongDNA{4}("ACGTACGT"), 3:5))
    "ACACGT"
```
"""
function delete_at!(seq :: LongDNA, pos_range :: UnitRange{Int})
    isempty(pos_range) && return seq
    return deleteat!(seq, pos_range)
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
        for _ in lastindex(seq):pos-1
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
    julia> String(copy_number_variation!(LongDNA{4}("ACGT"), 2, (LongDNA{4}("GG"), 3)))
    "AGGGGGGCGT"
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
       |- pos, deltatype, index
       |- pos, deltatype, index
```
"""
function _print_node(node::JSTNode, indent::Int)
    # build a label from its delta dict
    first_e  = first(values(node.delta))
    pos      = first_e.position
    idxs     = sort!(collect(keys(node.delta)))
    dtype    = first_e.delta_type
    label    = "pos=$pos idxs=$(idxs) type=$(dtype)"

    # print with the right indent
    println(" "^ (2*indent) * "|- " * label)

    # recurse into children
    for child in node.children
        _print_node(child, indent+1)
    end
end

# top-level: print each root child
function print_tree(tree :: JSTree)
    build_tree!(tree)
    for node in tree.rootChildren
        _print_node(node, 0)
    end
    return nothing
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
    io = IOBuffer()
    for i in 1:length(js.deltaMap)
        modified_seq = apply_delta(js.reference, js.deltaMap[i])
        println(io, "Sequence ", i, ": ", modified_seq)
    end
    return String(take!(io))
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
function _check_indices(n :: Integer, indices)
    for idx in indices
        1 <= idx <= n || throw(BoundsError(1:n, idx))
    end
    return nothing
end

function add_delta!(js :: JournaledString,
                    indices :: Vector{Int},
                    delta_type :: DeltaType,
                    position :: Int,
                    data)
    _check_indices(length(js.deltaMap), indices)
    for idx in indices
        new_entry = JournalEntry(delta_type, position, data, js.current_time)
        js.deltaMap[idx][js.current_time] = new_entry
        js.current_time += UInt64(1)
    end
    return nothing
end

"""
Return the bucket of entries that share `delta_type` and `data` at one
position, creating it when no such bucket exists.

Buckets are matched by value, so two separately built but equal payloads
land in the same branch.
"""
function _matching_bucket!(buckets :: Vector{Dict{Int64, JournalEntry}},
                           delta_type :: DeltaType,
                           data)
    for bucket in buckets
        isempty(bucket) && continue
        entry = first(values(bucket))
        if entry.delta_type == delta_type && entry.payload == make_payload(delta_type, data)
            return bucket
        end
    end
    bucket = Dict{Int64, JournalEntry}()
    push!(buckets, bucket)
    return bucket
end

"""
Adds a delta to every listed sequence in a JSTree.

Entries that share a position, delta type, and payload share one bucket,
which `build_tree!` later turns into one branch.
"""
function _insert_tree_entry!(tree :: JSTree, idx :: Integer, entry :: JournalEntry)
    vec = get!(tree.journal, entry.position) do
        Vector{Dict{Int64, JournalEntry}}()
    end
    bucket = _matching_bucket!(vec, entry.delta_type, entry.payload)
    bucket[Int64(idx)] = entry
    return nothing
end

function add_delta!(tree :: JSTree,
                    indices :: Vector{Int},
                    delta_type :: DeltaType,
                    position :: Int,
                    data)
    _check_indices(tree.length, indices)
    for idx in indices
        entry = JournalEntry(delta_type, position, data, tree.current_time)
        _insert_tree_entry!(tree, idx, entry)
        tree.current_time += UInt64(1)
    end
    build_tree!(tree)
    return nothing
end

"""
Insert a JournalEntry template into a JSTree, stamping each insertion
with the tree's current_time (and then incrementing it).

# Args
- `tree`     : the JSTree to modify
- `indices`  : Vector of sequence indices where this entry applies
- `entry`    : a JournalEntry whose `delta_type`, `position`, and `data`
               we'll reuse (its `.time` is ignored)

# Returns
- `nothing`  : mutates `tree.journal` and `tree.current_time`
"""
function add_delta!(tree :: JSTree, indices :: Vector{Int}, entry :: JournalEntry)
    _check_indices(tree.length, indices)
    for idx in indices
        stamped = JournalEntry(entry.delta_type, entry.position, entry.payload, tree.current_time)
        _insert_tree_entry!(tree, idx, stamped)
        tree.current_time += UInt64(1)
    end
    build_tree!(tree)
    return nothing
end

function build_tree!(tree::JSTree)
    tree.rootChildren = JSTNode[]             # clear any prior structure
    lastPrimary = nothing                       # tracks the main chain

    for (pos, buckets) in tree.journal
        if isempty(buckets)
            continue
        end

        # first bucket = primary chain node
        b0 = buckets[1]
        node0 = JSTNode(nothing, copy(b0), JSTNode[])
        if lastPrimary === nothing
            push!(tree.rootChildren, node0)
        else
            node0.parent = lastPrimary
            push!(lastPrimary.children, node0)
        end
        lastPrimary = node0

        # any other buckets at same pos = secondary branches
        for b in buckets[2:end]
            node = JSTNode(nothing, copy(b), JSTNode[])
            parentNode = node0.parent
            node.parent = parentNode
            if parentNode === nothing
                push!(tree.rootChildren, node)
            else
                push!(parentNode.children, node)
            end
        end
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
    _check_indices(length(js.deltaMap), indices)
    for idx in indices
        stamped = JournalEntry(entry.delta_type, entry.position, entry.payload, js.current_time)
        js.deltaMap[idx][js.current_time] = stamped
        js.current_time += UInt64(1)
    end
    return nothing
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
function remove_delta!(js :: JournaledString, idx :: Int, time :: Integer)
    _check_indices(length(js.deltaMap), (idx,))
    key = Timestamp(time)
    delta = js.deltaMap[idx]
    if !haskey(delta, key)
        throw(ArgumentError("No mutation found at time: $time at index: $idx"))
    end
    delete!(delta, key)
    return nothing
end

"""
    remove_delta!(tree::JSTree, time)

Remove the journal entry stamped with `time` and rebuild the derived tree.
"""
function remove_delta!(tree::JSTree, time::Timestamp)
    found = false
    pos_to_clean = nothing
    bucket_idx_to_clean = nothing
    
    # Search through all journal entries
    for (pos, buckets) in tree.journal
        for (bucket_idx, bucket) in enumerate(buckets)
            for (idx, entry) in bucket
                if entry.time == time
                    # Found target entry - remove from bucket
                    delete!(bucket, idx)
                    found = true
                    pos_to_clean = pos
                    bucket_idx_to_clean = bucket_idx
                    break
                end
            end
            found && break
        end
        found && break
    end

    if !found
        throw(ArgumentError("No mutation found at time: $time"))
    end

    # Clean up empty structures
    if !isnothing(pos_to_clean)
        bucket = tree.journal[pos_to_clean][bucket_idx_to_clean]
        
        # Remove bucket if empty
        if isempty(bucket)
            deleteat!(tree.journal[pos_to_clean], bucket_idx_to_clean)
            
            # Remove position entry if no buckets left
            if isempty(tree.journal[pos_to_clean])
                delete!(tree.journal, pos_to_clean)
            end
        end
    end
    build_tree!(tree)
    return nothing
end

remove_delta!(tree :: JSTree, time :: Integer) = remove_delta!(tree, Timestamp(time))

"""
Count total deltas in journal (for testing)
"""
function delta_count(tree::JSTree)
    count = 0
    for (_, buckets) in tree.journal
        for bucket in buckets
            count += length(bucket)
        end
    end
    return count
end


### functions.jl ends here.
