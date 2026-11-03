### -*- Mode: Julia -*-

### functions.jl


"""
    delete_at!(seq, range) -> LongSequence

Delete `range` from `seq` and return `seq`. An empty range leaves `seq`
unchanged. The range is 1-based and closed, as in `deleteat!`.
"""
function delete_at!(seq :: LongSequence, pos_range :: UnitRange{<:Integer})
    isempty(pos_range) && return seq
    return deleteat!(seq, pos_range)
end


function _append_gaps!(seq :: LongSequence, n :: Integer)
    n <= 0 && return seq
    symbol = gap(eltype(seq))
    symbol in symbols(Alphabet(seq)) ||
        throw(ArgumentError(
            "alphabet $(typeof(Alphabet(seq))) cannot encode a gap; " *
            "use a 4-bit nucleic alphabet such as LongDNA{4} or LongRNA{4} " *
            "for a structural variant past the reference end",
        ))
    for _ in 1:n
        push!(seq, symbol)
    end
    return seq
end

"""
    structure_variation!(seq, pos, allele) -> LongSequence

Splice `allele` into `seq` at `pos` and return `seq`.

A position before the last base inserts `allele` immediately before that
base. The last position appends `allele`. A position past the end inserts
`pos - lastindex(seq)` gap symbols and then appends `allele`. `pos` must
be at least 1. A 2-bit alphabet cannot encode the gap pad.
"""
function structure_variation!(seq :: LongSequence, pos :: Integer, subseq :: LongSequence)
    pos >= 1 || throw(ArgumentError("position must be >= 1, got $pos"))
    if pos == lastindex(seq)
        append!(seq, subseq)
    elseif pos < lastindex(seq)
        spliceinto!(seq, pos, subseq)
    elseif pos > lastindex(seq)
        _append_gaps!(seq, pos - lastindex(seq))
        append!(seq, subseq)
    end
    return seq
end


"""
    copy_number_variation!(seq, pos, (allele, copies)) -> LongSequence

Insert `copies` tandem copies of `allele` at `pos` and return `seq`.

Each copy is spliced at the same position, so identical copies stack to
the left. `pos` and `copies` must be at least 1.
"""
function copy_number_variation!(seq :: LongSequence,
                                pos :: Integer,
                                params :: Tuple{<:LongSequence, <:Integer})
    subseq, copies = params
    pos >= 1 || throw(ArgumentError("position must be >= 1, got $pos"))
    copies >= 1 || throw(ArgumentError("copy number must be positive, got $copies"))
    for _ in 1:copies
        spliceinto!(seq, pos, subseq)
    end
    return seq
end


### End of deltas

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

"""
    print_tree(tree)

Print the tree derived from the journal of `tree`.

Each line is one edit bucket: reference position, sequence indices, and
delta type. The tree is rebuilt from the journal before printing. Output
goes to `stdout`.
"""
function print_tree(tree :: JSTree)
    build_tree!(tree)
    for node in tree.rootChildren
        _print_node(node, 0)
    end
    return nothing
end

"""
    print_sequences(journal)

Print every derived sequence of `journal` to `stdout`, one per line,
as `Sequence <index>: <sequence>`.
"""
function print_sequences(jst::JournaledString)
    for i in 1:length(jst.deltaMap)
        modified_seq = apply_delta(jst.reference, jst.deltaMap[i])
        println("Sequence $i: ", modified_seq)
    end
end


"""
    build_sequences(journal) -> String

Return the text printed by [`print_sequences`](@ref), without writing
to `stdout`.
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
    print_deltas(journal)

Print every edit in `journal` to `stdout`, grouped by sequence index.
Each edit line is `[time = <timestamp>] JournalEntry: <entry>`.
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


function _check_indices(n :: Integer, indices)
    for idx in indices
        1 <= idx <= n || throw(BoundsError(1:n, idx))
    end
    return nothing
end

"""
    add_delta!(journal, indices, kind, position, data)
    add_delta!(journal, indices, entry)
    add_delta!(journal, index, kind, position, data)
    add_delta!(journal, index, entry)

Record one edit on each listed derived sequence of a [`JournaledString`](@ref)
or [`JSTree`](@ref).

`indices` is a 1-based vector of sequence indexes. A bare integer records
the edit on that one sequence. An index outside the journal raises
`BoundsError`.

`position` is a 1-based reference coordinate. `kind` is a [`DeltaType`](@ref)
and `data` is the payload accepted by [`JournalEntry`](@ref). Passing an
`entry` reuses its kind, position, and payload. The timestamp already stored
on `entry` is ignored.

Each recorded edit is stamped with the journal's current clock, and the
clock then advances by one. A call that lists several sequences stores the
same allele at successive times. The allele is checked against the reference
before it is stored. A deletion past the end of the reference, a symbol
outside the alphabet, or a structural variant that needs a gap on a 2-bit
alphabet raises `ArgumentError`.

# Examples
```jldoctest
julia> journal = JournaledString(dna"ACGTACGT", 1);

julia> add_delta!(journal, 1, DeltaTypeDel, 3, 2);

julia> String(sequence(journal, 1))
"ACACGT"
```
"""
function add_delta!(js :: JournaledString{S},
                    indices :: AbstractVector{<:Integer},
                    delta_type :: DeltaType,
                    position :: Integer,
                    data) where {S <: LongSequence}
    _check_indices(length(js.deltaMap), indices)
    for idx in indices
        new_entry = JournalEntry{S}(delta_type, position, data, js.current_time)
        _ensure_applicable(js.reference, new_entry)
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
function _matching_bucket!(::Type{S},
                           buckets :: Vector{Dict{Int64, JournalEntry{S}}},
                           delta_type :: DeltaType,
                           data) where {S <: LongSequence}
    payload = make_payload(S, delta_type, data)
    for bucket in buckets
        isempty(bucket) && continue
        entry = first(values(bucket))
        if entry.delta_type == delta_type && entry.payload == payload
            return bucket
        end
    end
    bucket = Dict{Int64, JournalEntry{S}}()
    push!(buckets, bucket)
    return bucket
end

"""
Place `entry` in the bucket of equal edits at its reference position.
"""
function _insert_tree_entry!(tree :: JSTree{S}, idx :: Integer, entry :: JournalEntry{S}) where {S <: LongSequence}
    vec = get!(tree.journal, entry.position) do
        Vector{Dict{Int64, JournalEntry{S}}}()
    end
    bucket = _matching_bucket!(S, vec, entry.delta_type, entry.payload)
    bucket[Int64(idx)] = entry
    return nothing
end

function add_delta!(tree :: JSTree{S},
                    indices :: AbstractVector{<:Integer},
                    delta_type :: DeltaType,
                    position :: Integer,
                    data) where {S <: LongSequence}
    _check_indices(tree.length, indices)
    for idx in indices
        entry = JournalEntry{S}(delta_type, position, data, tree.current_time)
        _ensure_applicable(tree.root, entry)
        _insert_tree_entry!(tree, idx, entry)
        tree.current_time += UInt64(1)
    end
    build_tree!(tree)
    return nothing
end

function add_delta!(tree :: JSTree{S}, indices :: AbstractVector{<:Integer}, entry :: JournalEntry) where {S <: LongSequence}
    _check_indices(tree.length, indices)
    for idx in indices
        stamped = JournalEntry{S}(entry.delta_type, entry.position, entry.payload, tree.current_time)
        _ensure_applicable(tree.root, stamped)
        _insert_tree_entry!(tree, idx, stamped)
        tree.current_time += UInt64(1)
    end
    build_tree!(tree)
    return nothing
end

function build_tree!(tree :: JSTree{S}) where {S <: LongSequence}
    tree.rootChildren = JSTNode{S}[]
    lastPrimary = nothing

    for (pos, buckets) in tree.journal
        if isempty(buckets)
            continue
        end

        # first bucket = primary chain node
        b0 = buckets[1]
        node0 = JSTNode{S}(nothing, copy(b0), JSTNode{S}[])
        if lastPrimary === nothing
            push!(tree.rootChildren, node0)
        else
            node0.parent = lastPrimary
            push!(lastPrimary.children, node0)
        end
        lastPrimary = node0

        # any other buckets at same pos = secondary branches
        for b in buckets[2:end]
            node = JSTNode{S}(nothing, copy(b), JSTNode{S}[])
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

function add_delta!(js :: JournaledString{S},
                    indices :: AbstractVector{<:Integer},
                    entry :: JournalEntry) where {S <: LongSequence}
    _check_indices(length(js.deltaMap), indices)
    for idx in indices
        stamped = JournalEntry{S}(entry.delta_type, entry.position, entry.payload, js.current_time)
        _ensure_applicable(js.reference, stamped)
        js.deltaMap[idx][js.current_time] = stamped
        js.current_time += UInt64(1)
    end
    return nothing
end

function add_delta!(journal :: Union{JournaledString, JSTree},
                    index :: Integer,
                    delta_type :: DeltaType,
                    position :: Integer,
                    data)
    return add_delta!(journal, [Int(index)], delta_type, position, data)
end

function add_delta!(journal :: Union{JournaledString, JSTree},
                    index :: Integer,
                    entry :: JournalEntry)
    return add_delta!(journal, [Int(index)], entry)
end

"""
    remove_delta!(journal, index, time)
    remove_delta!(tree, time)

Remove one recorded edit.

On a [`JournaledString`](@ref), `index` selects the derived sequence and
`time` selects the edit. On a [`JSTree`](@ref), `time` selects the edit
and the derived tree is rebuilt. A missing edit raises `ArgumentError`.
An index outside the journal raises `BoundsError`.
"""
function remove_delta!(js :: JournaledString, idx :: Integer, time :: Integer)
    _check_indices(length(js.deltaMap), (idx,))
    key = Timestamp(time)
    delta = js.deltaMap[idx]
    if !haskey(delta, key)
        throw(ArgumentError("No mutation found at time: $time at index: $idx"))
    end
    delete!(delta, key)
    return nothing
end

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
    delta_count(tree) -> Int

Number of edits stored in `tree`.
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
