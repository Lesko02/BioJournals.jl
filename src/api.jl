### -*- Mode: Julia -*-

### api.jl


"""
    sequence(journal, i) -> LongSequence

Materialize derived sequence `i` by replaying its edits on the reference.

This method belongs to `FASTX.sequence`, so a journal and a FASTA or FASTQ
record share one name. The result is an ordinary `LongSequence`. BioSequences
search and translation, and the FASTX writers, accept it directly.

`i` is 1-based. A [`JournaledString`](@ref) indexes `deltaMap`. A
[`JSTree`](@ref) indexes the sequences tracked by the tree. An index outside
that range raises `BoundsError`.

# Examples
```jldoctest
julia> journal = JournaledString(dna"ACGT", 1);

julia> add_delta!(journal, 1, DeltaTypeIns, 2, "TT");

julia> sequence(journal, 1)
6nt DNA Sequence:
ATTCGT
```
"""
function sequence(js :: JournaledString, i :: Integer)
    return apply_delta(js.reference, edits_at(js, i))
end

function sequence(tree :: JSTree, i :: Integer)
    return apply_delta(tree.root, edits_at(tree, i))
end

"""
    sequences(journal) -> Vector{<:LongSequence}

Materialize every derived sequence in `journal`, in index order.
"""
function sequences(js :: JournaledString)
    return [sequence(js, i) for i in eachindex(js.deltaMap)]
end

function sequences(tree :: JSTree)
    return [sequence(tree, i) for i in 1:tree.length]
end

"""
    nsequences(journal) -> Int

Number of derived sequences stored in `journal`.
"""
nsequences(js :: JournaledString) = length(js.deltaMap)
nsequences(tree :: JSTree) = tree.length

"""
    get_mutation_history(delta_map) -> String

List every edit in `delta_map`, oldest first. Each line has the form
`Time <timestamp>: <entry>`.
"""
function get_mutation_history(delta_map :: DeltaMap)
    io = IOBuffer()
    for (time, entry) in delta_map
        println(io, "Time ", UInt64(time), ": ", entry)
    end
    return String(take!(io))
end


"""
    get_mutation_interval(journal, time1, time2) -> String

List edits whose timestamp lies in the closed interval `time1:time2`.

Each line names the time, the entry, and the 1-based sequence index.
An empty string means no edit falls in the interval.
"""
function get_mutation_interval(js :: JournaledString,
                               time1 :: Int,
                               time2 :: Int)
    
    io = IOBuffer()
    for i in 1:length(js.deltaMap)
        for (time, entry) in js.deltaMap[i]
            if time1 <= time <= time2
                println(io, "Time ", time, ": ", entry, " at index ", i)
            end
        end
    end
    return String(take!(io))
end


"""
    get_sequences_at_time(journal, time) -> Vector{<:LongSequence}

Materialize every derived sequence using only edits whose timestamp is
less than or equal to `time`. Later edits stay in the journal and are
omitted from the result. A negative `time` raises `DomainError`.
"""
function get_sequences_at_time(jst :: JournaledString{S}, time :: Integer) where {S <: LongSequence}
    limit = Timestamp(time)
    sequences_at_time = Vector{S}(undef, length(jst.deltaMap))
    for i in eachindex(jst.deltaMap)
        kept = JournalEntry{S}[]
        for (entry_time, entry) in jst.deltaMap[i]
            entry_time > limit && break
            push!(kept, entry)
        end
        sequences_at_time[i] = apply_delta(jst.reference, kept)
    end
    return sequences_at_time
end


"""
    simulate_mutation!(journal, index, entry)

Record `entry` on derived sequence `index`.

The stored edit receives `journal.current_time`. The timestamp on `entry`
is ignored. This is [`add_delta!`](@ref) for one sequence.
"""
function simulate_mutation!(jst :: JournaledString,
                            index :: Integer,
                            entry :: JournalEntry)
   add_delta!(jst, [index], entry)
end


"""
    remove_mutation!(journal, index, time)

Remove the edit stored on sequence `index` at `time`.

Raises `ArgumentError` when that timestamp is absent, and `BoundsError`
when `index` is outside the journal.
"""
function remove_mutation!(jst :: JournaledString,
                          index :: Integer,
                          time :: Integer)
    remove_delta!(jst ,index, time)
end


"""
    is_equal(a, b) -> Bool

Compare the contents of two journals.

Two [`JournaledString`](@ref) values are equal when their reference,
clock, reference identifier, and per-sequence journals match. Two
[`JSTree`](@ref) values are equal when their reference, clock, number of
sequences, and journal match. The derived `rootChildren` cache is ignored.
"""
function is_equal(jst1 :: JournaledString,
                  jst2 :: JournaledString) :: Bool
    return jst1.reference == jst2.reference &&
           jst1.current_time == jst2.current_time &&
           jst1.reference_id == jst2.reference_id &&
           jst1.deltaMap == jst2.deltaMap
end


function is_equal(jst1 :: JSTree, jst2 :: JSTree) :: Bool
    return jst1.root == jst2.root &&
           jst1.current_time == jst2.current_time &&
           jst1.length == jst2.length &&
           jst1.journal == jst2.journal
end



"""
    print_results(results)

Print the matches returned by [`exact_search`](@ref) or
[`approximate_search`](@ref).

Integer keys are sequence indexes. String keys are labels. Keys are printed
in sorted order. An empty vector is reported as no match. The function
writes to `stdout` and returns `nothing`.
"""
function print_results(results :: Dict{Int64, Vector{UnitRange{Int64}}})
    for i in sort!(collect(keys(results)))
        ranges = results[i]
        if isempty(ranges)
            println("No Match in DeltaMap $i")
        else
            println("Match in $i:")
            println("Ranges: ", ranges)
        end
    end
    return nothing
end


function print_results(results :: Dict{String, Vector{UnitRange{Int64}}})
    for name in sort!(collect(keys(results)))
        ranges = results[name]
        if isempty(ranges)
            println("No Match at $name")
        else
            println("Match at $name:")
            println("Ranges: ", ranges)
        end
    end
    return nothing
end


### api.jl ends here.
