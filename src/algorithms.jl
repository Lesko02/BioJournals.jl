### -*- Mode: Julia -*-

### algorithms.jl
#
# Exact and approximate search share one path. The reference is scanned
# once. Sequences with no edits reuse those hits. Sequences with edits keep
# the hits that no edit touches, mapped into the edited coordinate system,
# and rescan only a window around each edit wide enough to contain every
# new or broken match.

function approximate_findall(query, tolerance :: Int, seq :: LongSequence)
    results = UnitRange{Int64}[]
    pos = findfirst(query, tolerance, seq)
    while pos !== nothing && !isempty(pos)
        push!(results, pos)
        pos = findnext(query, tolerance, seq, last(pos) + tolerance + 1)
    end
    return results
end

function _check_tolerance(tol :: Integer)
    1 <= tol <= 99 ||
        throw(ArgumentError("tolerance must be between 1 and 99 percent, got $tol"))
    return nothing
end

_error_budget(needle, tol :: Integer) = ceil(Int, length(needle) * tol / 100)

function _edit_touches_span(entry :: JournalEntry, span :: UnitRange{Int}, ref_len :: Int)
    position = entry.position
    payload = entry.payload
    if payload isa Substitution
        return first(span) <= position <= last(span)
    elseif payload isa Deletion
        return position <= last(span) && position + payload.len - 1 >= first(span)
    elseif payload isa StructuralVariant && position >= ref_len
        return last(span) == ref_len
    else
        # An insertion anchored at the first base of the span sits in front of
        # it and only shifts the hit. An anchor inside the span rewrites it.
        return first(span) < position <= last(span)
    end
end

function _span_affected(span :: UnitRange{Int}, edits, ref_len :: Int)
    for entry in edits
        _edit_touches_span(entry, span, ref_len) && return true
    end
    return false
end

"""Map a reference interval that no edit rewrites into the edited sequence."""
function _project_span(span :: UnitRange{Int}, edits, ref_len :: Int)
    start = _ref_base_pos(edits, first(span), ref_len)
    stop = _ref_base_pos(edits, last(span), ref_len)
    return Int64(start):Int64(stop)
end

function _occupied_range(entry :: JournalEntry, edits, ref_len :: Int, seq_len :: Int)
    payload = entry.payload
    if payload isa StructuralVariant && entry.position >= ref_len
        # Include the last reference base so a hit that ends there is
        # inside the rescan window, not only the appended tail.
        anchor = max(1, seq_len - _terminal_tail(edits, ref_len))
        return anchor:seq_len
    end
    if _anchor_deleted(edits, entry.position) && !(payload isa Deletion)
        return 1:0
    end
    if payload isa Deletion
        next_ref = entry.position + payload.len
        junction = if next_ref <= ref_len
            _ref_base_pos(edits, next_ref, ref_len)
        elseif entry.position == 1
            1
        else
            _ref_base_pos(edits, entry.position - 1, ref_len) + 1
        end
        return junction:junction
    elseif payload isa Substitution
        base = _ref_base_pos(edits, entry.position, ref_len)
        return base:base
    else
        base = _ref_base_pos(edits, entry.position, ref_len)
        lead = _anchor_insert_len(edits, entry.position, ref_len)
        start = base - lead
        stop = base - 1
        return start <= stop ? (start:stop) : (base:base)
    end
end

function _terminal_tail(edits, ref_len :: Int)
    total = 0
    for entry in edits
        payload = entry.payload
        payload isa StructuralVariant || continue
        entry.position >= ref_len || continue
        total += length(payload.seq) + (entry.position - ref_len)
    end
    return total
end

function _search_window(occupied :: UnitRange{Int}, context :: Int, seq_len :: Int)
    isempty(occupied) && return 1:0
    lo = max(1, first(occupied) - context + 1)
    hi = min(seq_len, last(occupied) + context - 1)
    return lo <= hi ? (lo:hi) : (1:0)
end

"""
Keep the same matches a left-to-right scanner would report.
`resume` is how many bases after the end of a hit the next search may start.
Exact search uses 1. Approximate search uses `tolerance + 1`.
"""
function _greedy_hits(candidates, resume :: Int)
    ordered = sort!(collect(Set{UnitRange{Int64}}(candidates)); by = r -> (first(r), last(r)))
    kept = UnitRange{Int64}[]
    next_ok = 1
    for hit in ordered
        first(hit) < next_ok && continue
        push!(kept, hit)
        next_ok = last(hit) + resume
    end
    return kept
end

function _hits_on_sequence(search_window! :: Function,
                           reference :: LongSequence,
                           edits,
                           all_ref_hits,
                           clean_hits,
                           context :: Int,
                           resume :: Int)
    isempty(edits) && return copy(clean_hits)
    ref_len = length(reference)
    seq = apply_delta(reference, edits)
    candidates = UnitRange{Int64}[]
    for span in all_ref_hits
        _span_affected(span, edits, ref_len) && continue
        projected = _project_span(span, edits, ref_len)
        length(projected) == length(span) && push!(candidates, projected)
    end
    seq_len = length(seq)
    seen = Set(candidates)
    for entry in edits
        occupied = _occupied_range(entry, edits, ref_len, seq_len)
        window = _search_window(occupied, context, seq_len)
        isempty(window) && continue
        search_window!(candidates, seen, seq, window)
    end
    return _greedy_hits(candidates, resume)
end

function _exact_window!(hits, seen, seq, window, query)
    offset = first(window) - 1
    for match in _all_exact(query, seq[window])
        hit = Int64(first(match) + offset):Int64(last(match) + offset)
        if hit ∉ seen
            push!(seen, hit)
            push!(hits, hit)
        end
    end
    return nothing
end

function _approx_window!(hits, seen, seq, window, query, budget)
    offset = first(window) - 1
    for match in _all_approx(query, budget, seq[window])
        hit = Int64(first(match) + offset):Int64(last(match) + offset)
        if hit ∉ seen
            push!(seen, hit)
            push!(hits, hit)
        end
    end
    return nothing
end

function _sequence_count(js :: JournaledString)
    return length(js.deltaMap)
end

function _sequence_count(tree :: JSTree)
    return tree.length
end

function _reference_of(js :: JournaledString)
    return js.reference
end

function _reference_of(tree :: JSTree)
    return tree.root
end

function _all_exact(query, seq :: LongSequence)
    hits = UnitRange{Int64}[]
    pos = findfirst(query, seq)
    while pos !== nothing && !isempty(pos)
        push!(hits, pos)
        nxt = first(pos) + 1
        nxt > lastindex(seq) && break
        pos = findnext(query, seq, nxt)
    end
    return hits
end

"""
    exact_search(journal, needle) -> Dict{Int, Vector{UnitRange{Int}}}

Find exact matches of `needle` in every derived sequence.

The result maps a 1-based sequence index to the match intervals in that
derived sequence. Intervals are 1-based and closed. A sequence with no
edits reuses the matches found on the reference.
"""
function exact_search(journal, needle :: LongSequence)
    reference = _reference_of(journal)
    query = ExactSearchQuery(needle)
    clean_hits = collect(BioSequences.findall(query, reference))
    all_hits = _all_exact(query, reference)
    n = _sequence_count(journal)
    context = max(length(needle), 1)
    results = Dict{Int64, Vector{UnitRange{Int64}}}()
    for index in 1:n
        edits = edits_at(journal, index)
        results[index] = _hits_on_sequence(
            reference, edits, all_hits, clean_hits, context, 1) do hits, seen, seq, window
            _exact_window!(hits, seen, seq, window, query)
        end
    end
    return results
end

"""
    approximate_search(journal, needle, tol=5) -> Dict{Int, Vector{UnitRange{Int}}}

Find `needle` in every sequence of a `JournaledString` or `JSTree`.

`tol` is the allowed mismatch percentage, from 1 through 99. The default
is 5. Hits are sorted and unique.
"""
function approximate_search(journal, needle :: LongSequence, tol :: Integer = 5)
    _check_tolerance(tol)
    reference = _reference_of(journal)
    budget = _error_budget(needle, tol)
    query = ApproximateSearchQuery(needle)
    clean_hits = approximate_findall(query, budget, reference)
    all_hits = _all_approx(query, budget, reference)
    n = _sequence_count(journal)
    # Indels inside an approximate hit can extend it by `budget` bases.
    context = max(length(needle) + 2 * budget, 1)
    resume = budget + 1
    results = Dict{Int64, Vector{UnitRange{Int64}}}()
    for index in 1:n
        edits = edits_at(journal, index)
        results[index] = _hits_on_sequence(
            reference, edits, all_hits, clean_hits, context, resume) do hits, seen, seq, window
            _approx_window!(hits, seen, seq, window, query, budget)
        end
    end
    return results
end

"""Every approximate hit, including ones the public scanner skips as overlaps."""
function _all_approx(query, budget :: Int, seq :: LongSequence)
    hits = UnitRange{Int64}[]
    pos = findfirst(query, budget, seq)
    while pos !== nothing && !isempty(pos)
        push!(hits, pos)
        nxt = first(pos) + 1
        nxt > lastindex(seq) && break
        pos = findnext(query, budget, seq, nxt)
    end
    return hits
end

### algorithms.jl ends here.
