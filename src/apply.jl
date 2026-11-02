### -*- Mode: Julia -*-

### apply.jl
#
# Journal entries store 1-based coordinates on the reference. Replay walks
# the edits in time order and maps each reference coordinate through the
# edits already applied, so deleting an earlier event does not require
# rebasing the ones that remain.

function _ins_len(entry :: JournalEntry, ref_len :: Int)
    payload = entry.payload
    if payload isa Insertion
        return length(payload.seq)
    elseif payload isa CopyNumber
        return length(payload.seq) * payload.copies
    elseif payload isa StructuralVariant && entry.position < ref_len
        return length(payload.seq)
    else
        return 0
    end
end

function _anchor_deleted(edits, anchor :: Int)
    for entry in edits
        payload = entry.payload
        payload isa Deletion || continue
        last_base = entry.position + payload.len - 1
        if entry.position <= anchor <= last_base
            return true
        end
    end
    return false
end

"""
How many sequence bases this edit inserts or deletes before reference base
`ref_index`. Insertions anchored on a base that a deletion has already
removed do not shift later coordinates.
"""
function _length_delta_before(entry :: JournalEntry, ref_index :: Int, ref_len :: Int, edits)
    payload = entry.payload
    if payload isa Deletion
        last_base = entry.position + payload.len - 1
        return ref_index > last_base ? -payload.len : 0
    end
    added = _ins_len(entry, ref_len)
    added == 0 && return 0
    entry.position <= ref_index || return 0
    return _anchor_deleted(edits, entry.position) ? 0 : added
end

function _ref_base_pos(edits, ref_index :: Int, ref_len :: Int)
    1 <= ref_index <= ref_len ||
        throw(ArgumentError("reference position $ref_index is outside 1:$ref_len"))
    if _anchor_deleted(edits, ref_index)
        throw(ArgumentError("reference position $ref_index falls inside a deletion"))
    end
    shift = 0
    for entry in edits
        shift += _length_delta_before(entry, ref_index, ref_len, edits)
    end
    return ref_index + shift
end

function _anchor_insert_len(edits, anchor :: Int, ref_len :: Int)
    _anchor_deleted(edits, anchor) && return 0
    total = 0
    for entry in edits
        entry.position == anchor || continue
        total += _ins_len(entry, ref_len)
    end
    return total
end

function _check_edit(entry :: JournalEntry, ref_len :: Int)
    payload = entry.payload
    if payload isa Deletion
        last_base = entry.position + payload.len - 1
        last_base <= ref_len ||
            throw(ArgumentError("deletion at $(entry.position) length $(payload.len) extends past reference length $ref_len"))
    elseif payload isa StructuralVariant
        return nothing
    elseif entry.position > ref_len
        throw(ArgumentError("position $(entry.position) is past reference length $ref_len"))
    end
    return nothing
end

function _repeat(seq :: LongDNA{4}, copies :: Int)
    copies == 1 && return copy(seq)
    block = LongDNA{4}()
    for _ in 1:copies
        append!(block, seq)
    end
    return block
end

"""
Apply `entry` to `seq`. `prior` contains the edits already applied, in time
order, and `ref_len` is the length of the untouched reference.
"""
function _apply_reference!(seq :: LongDNA{4}, entry :: JournalEntry, prior, ref_len :: Int)
    _check_edit(entry, ref_len)
    payload = entry.payload
    if payload isa Insertion
        insert!(seq, _ref_base_pos(prior, entry.position, ref_len), payload.seq)
    elseif payload isa Substitution
        seq[_ref_base_pos(prior, entry.position, ref_len)] = payload.base
    elseif payload isa CopyNumber
        block = _repeat(payload.seq, payload.copies)
        insert!(seq, _ref_base_pos(prior, entry.position, ref_len), block)
    elseif payload isa Deletion
        last_ref = entry.position + payload.len - 1
        start_base = _ref_base_pos(prior, entry.position, ref_len)
        lead = _anchor_insert_len(prior, entry.position, ref_len)
        end_base = _ref_base_pos(prior, last_ref, ref_len)
        delete_at!(seq, (start_base - lead):end_base)
    elseif payload isa StructuralVariant
        if entry.position >= ref_len
            gaps = entry.position - ref_len
            if gaps > 0
                append!(seq, LongDNA{4}("-"^gaps))
            end
            append!(seq, payload.seq)
        else
            insert!(seq, _ref_base_pos(prior, entry.position, ref_len), payload.seq)
        end
    else
        throw(ArgumentError("unsupported payload $(typeof(payload))"))
    end
    return seq
end

"""
    apply_delta(reference, entries) -> LongDNA{4}

Replay `entries` against `reference`. Each position is a reference
coordinate. Edits are applied in timestamp order.
"""
function apply_delta(reference :: LongDNA{4}, entries :: AbstractVector{JournalEntry})
    seq = copy(reference)
    isempty(entries) && return seq
    ordered = sort(entries; by = entry -> UInt64(entry.time))
    ref_len = length(reference)
    prior = JournalEntry[]
    sizehint!(prior, length(ordered))
    for entry in ordered
        _apply_reference!(seq, entry, prior, ref_len)
        push!(prior, entry)
    end
    return seq
end

apply_delta(reference :: LongDNA{4}, entry :: JournalEntry) =
    apply_delta(reference, JournalEntry[entry])

apply_delta(reference :: LongDNA{4}, delta :: DeltaMap) =
    apply_delta(reference, collect(values(delta)))

"""Edits for sequence `index`, oldest first."""
function edits_at(js :: JournaledString, index :: Integer)
    _check_indices(length(js.deltaMap), (index,))
    return collect(values(js.deltaMap[index]))
end

function edits_at(tree :: JSTree, index :: Integer)
    _check_indices(tree.length, (index,))
    found = JournalEntry[]
    for (_, buckets) in tree.journal
        for bucket in buckets
            if haskey(bucket, Int64(index))
                push!(found, bucket[Int64(index)])
            end
        end
    end
    sort!(found; by = entry -> UInt64(entry.time))
    return found
end

### apply.jl ends here.
