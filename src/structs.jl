### -*- Mode: Julia -*-

### structs.jl

"""
    Timestamp(x::Integer) -> Timestamp

Clock value stored on a [`JournalEntry`](@ref). The value is an unsigned
64-bit integer.

A negative `x` raises `DomainError`. Timestamps compare with other
timestamps and with integers, including `BigInt`. Adding an unsigned
integer advances the clock.

# Examples
```jldoctest
julia> Timestamp(3) < Timestamp(4)
true

julia> Timestamp(3) == 3
true
```
"""
primitive type Timestamp <: Unsigned 64 end

# `reinterpret` requires identical width, so every other integer converts
# through `UInt64` first.
Timestamp(x::UInt64) = reinterpret(Timestamp, x)
function Timestamp(x::Integer)
    x < 0 && throw(DomainError(x, "Timestamp must be non-negative"))
    return Timestamp(UInt64(x))
end

Base.UInt64(ts::Timestamp) = reinterpret(UInt64, ts)
Base.convert(::Type{Timestamp}, x::Integer) = Timestamp(x)
Base.convert(::Type{Int64}, ts::Timestamp) = reinterpret(Int64, ts)
addto(ts::Timestamp, x::Unsigned) = Timestamp(UInt64(ts) + x)
addto(ts::Timestamp, x::Timestamp) = Timestamp(UInt64(ts) + UInt64(x))

Base.:+(ts::Timestamp, n::Unsigned) = Timestamp(UInt64(ts) + n)
Base.:+(ts::Timestamp, n::UInt64) = Timestamp(UInt64(ts) + n)

Base.:(==)(a::Timestamp, b::Timestamp) = UInt64(a) == UInt64(b)
Base.isequal(a::Timestamp, b::Timestamp) = UInt64(a) == UInt64(b)
Base.hash(ts::Timestamp, h::UInt64) = hash(UInt64(ts), h)
Base.isless(a::Timestamp, b::Timestamp) = UInt64(a) < UInt64(b)

Base.:(<)(a::Timestamp, b::Timestamp) = isless(a, b)
Base.:(<=)(a::Timestamp, b::Timestamp) = UInt64(a) <= UInt64(b)
Base.:(==)(ts::Timestamp, x::Integer) = x >= 0 && UInt64(ts) == UInt64(x)
Base.:(==)(x::Integer, ts::Timestamp) = ts == x
Base.:(<)(ts::Timestamp, x::Integer) = x >= 0 && UInt64(ts) < UInt64(x)
Base.:(<)(x::Integer, ts::Timestamp) = x < 0 || UInt64(x) < UInt64(ts)
Base.:<=(ts::Timestamp, x::Integer) = x >= 0 && UInt64(ts) <= UInt64(x)
Base.:<=(x::Integer, ts::Timestamp) = x < 0 || UInt64(x) <= UInt64(ts)

# `Timestamp <: Integer`, so the methods above are ambiguous with Base's
# `BigInt` methods. These overloads are the more specific pair.
Base.:(==)(ts::Timestamp, x::BigInt) = BigInt(UInt64(ts)) == x
Base.:(==)(x::BigInt, ts::Timestamp) = ts == x
Base.:(<)(ts::Timestamp, x::BigInt) = BigInt(UInt64(ts)) < x
Base.:(<)(x::BigInt, ts::Timestamp) = x < BigInt(UInt64(ts))
Base.:<=(ts::Timestamp, x::BigInt) = BigInt(UInt64(ts)) <= x
Base.:<=(x::BigInt, ts::Timestamp) = x <= BigInt(UInt64(ts))

Base.show(io::IO, ts::Timestamp) = show(io, UInt64(ts))
Base.print(io::IO, ts::Timestamp) = print(io, repr(UInt64(ts)))


"""
    DeltaType

Kind of edit stored in a [`JournalEntry`](@ref).

| Value | Payload | Effect |
|:--|:--|:--|
| `DeltaTypeDel` | positive length, or a [`Deletion`](@ref) | remove reference bases |
| `DeltaTypeIns` | sequence, or an [`Insertion`](@ref) | insert before a reference base |
| `DeltaTypeSnp` | one symbol, or a [`Substitution`](@ref) | replace one reference base |
| `DeltaTypeSV` | sequence, or a [`StructuralVariant`](@ref) | insert before a base, or append past the reference end |
| `DeltaTypeCNV` | `(sequence, copies)`, or a [`CopyNumber`](@ref) | insert tandem copies before a reference base |
"""
@enum DeltaType DeltaTypeDel DeltaTypeIns DeltaTypeSnp DeltaTypeSV DeltaTypeCNV

"""
    DeltaTypeDel

Remove one or more reference bases, starting at the entry position.
"""
DeltaTypeDel

"""
    DeltaTypeIns

Insert a sequence immediately before a reference base.
"""
DeltaTypeIns

"""
    DeltaTypeSnp

Replace the reference symbol at the entry position.
"""
DeltaTypeSnp

"""
    DeltaTypeSV

Insert a structural-variant allele before a reference base, or append it
when the position is at or beyond the last reference base.
"""
DeltaTypeSV

"""
    DeltaTypeCNV

Insert tandem copies of a sequence immediately before a reference base.
"""
DeltaTypeCNV

"""
    ReferenceInterval(first, last)

1-based closed interval on the reference. Both `first` and `last` are
included, matching `GenomicFeatures.Interval`.

`length` is `last - first + 1`. A position `p` is `in` the interval when
`first <= p <= last`.

An empty interval (`last == first - 1`) anchors an insertion immediately
before reference base `first`. The VCF record for that insertion has
`POS == first - 1`. `first` must be at least 1, and `last` may be at most
one base before `first`.

# Examples
```jldoctest
julia> length(ReferenceInterval(2, 4))
3

julia> isempty(ReferenceInterval(2, 1))
true

julia> 3 in ReferenceInterval(2, 4)
true
```
"""
struct ReferenceInterval
    first :: Int
    last :: Int

    function ReferenceInterval(first :: Integer, last :: Integer)
        start = Int(first)
        stop = Int(last)
        start >= 1 || throw(ArgumentError("interval start must be >= 1, got $start"))
        stop >= start - 1 ||
            throw(ArgumentError("interval end $stop is more than one base before start $start"))
        return new(start, stop)
    end
end

Base.first(interval :: ReferenceInterval) = interval.first
Base.last(interval :: ReferenceInterval) = interval.last
Base.isempty(interval :: ReferenceInterval) = interval.last < interval.first
Base.length(interval :: ReferenceInterval) = isempty(interval) ? 0 : interval.last - interval.first + 1
Base.in(position :: Integer, interval :: ReferenceInterval) = interval.first <= position <= interval.last
Base.:(==)(a :: ReferenceInterval, b :: ReferenceInterval) = a.first == b.first && a.last == b.last
Base.hash(interval :: ReferenceInterval, h :: UInt) = hash(interval.last, hash(interval.first, h))

function Base.show(io :: IO, interval :: ReferenceInterval)
    print(io, isempty(interval) ? "ReferenceInterval($(interval.first), before)" :
          "ReferenceInterval($(interval.first), $(interval.last))")
end

"""
    Deletion(len)

Payload that removes `len` reference bases, starting at the entry position.

`len` is a positive integer. The covered locus is
`ReferenceInterval(position, position + len - 1)`.
"""
struct Deletion
    len :: Int
end

"""
    Insertion(seq)

Payload that inserts `seq` immediately before the reference base named by
the entry position. `seq` uses the journal alphabet and is copied into the
journal. An empty sequence is rejected.
"""
struct Insertion{S <: LongSequence}
    seq :: S
end

"""
    Substitution(base)

Payload that replaces one reference symbol with `base`. `base` must belong
to the journal alphabet.
"""
struct Substitution{T}
    base :: T
end

"""
    StructuralVariant(seq)

Payload of a structural-variant edit. `seq` is copied into the journal and
must be non-empty.

When the entry position names an existing reference base, `seq` is inserted
immediately before that base. A position `p` at or beyond the last reference
base appends `seq`. The number of leading gap symbols is
`p - length(reference)`. A positive pad requires `LongDNA{4}` or
`LongRNA{4}`. `LongDNA{2}` and `LongRNA{2}` have no gap symbol, and
recording that edit raises `ArgumentError`.
"""
struct StructuralVariant{S <: LongSequence}
    seq :: S
end

"""
    CopyNumber(seq, copies)

Payload that inserts `copies` tandem copies of `seq` immediately before
the reference base named by the entry position. `copies` is a positive
integer and `seq` must be non-empty. Both are stored by copy.
"""
struct CopyNumber{S <: LongSequence}
    seq :: S
    copies :: Int
end

"""
    DeltaPayload{S}

Edit stored beside a reference of type `S`: a [`Deletion`](@ref),
[`Insertion`](@ref), [`Substitution`](@ref), [`StructuralVariant`](@ref),
or [`CopyNumber`](@ref).
"""
const DeltaPayload{S <: LongSequence} = Union{
    Deletion,
    Insertion{S},
    Substitution,
    StructuralVariant{S},
    CopyNumber{S},
}

function _as_seq(::Type{S}, seq :: S) where {S <: LongSequence}
    return copy(seq)
end
function _as_seq(::Type{S}, seq :: AbstractString) where {S <: LongSequence}
    return S(seq)
end
function _as_seq(::Type{S}, seq :: BioSequence) where {S <: LongSequence}
    return S(seq)
end

"""Element type of a long-sequence type. `eltype` on the type itself is `Any`."""
_symbol_type(::Type{LongSequence{A}}) where {A} = eltype(A)

function _as_symbol(::Type{S}, base) where {S <: LongSequence}
    symbol_type = _symbol_type(S)
    return base isa symbol_type ? base : convert(symbol_type, base)
end

function _sequence_type_of(data)
    if data isa Union{Insertion, StructuralVariant, CopyNumber}
        return typeof(data.seq)
    elseif data isa LongSequence
        return typeof(data)
    else
        return LongDNA{4}
    end
end

"""
    make_payload(kind, data)
    make_payload(::Type{S}, kind, data)

Translate a constructor argument into a [`DeltaPayload`](@ref) for a
reference of type `S`. The two-argument method uses `LongDNA{4}`.

Accepted inputs, by delta type:

  - deletion: a positive length, or a `Deletion`
  - insertion: a sequence, or an `Insertion`
  - SNP: a symbol of the reference alphabet, a `Char`, or a `Substitution`
  - structural variant: a sequence, or a `StructuralVariant`
  - copy number: `(sequence, copies)`, or a `CopyNumber`
"""
make_payload(kind :: DeltaType, data) = make_payload(LongDNA{4}, kind, data)

function make_payload(::Type{S}, kind :: DeltaType, data) where {S <: LongSequence}
    if kind == DeltaTypeDel
        data isa Deletion && return data
        len = Int(data)
        len >= 1 || throw(ArgumentError("deletion length must be positive, got $len"))
        return Deletion(len)
    elseif kind == DeltaTypeIns
        if data isa Insertion{S}
            isempty(data.seq) && throw(ArgumentError("insertion sequence must not be empty"))
            return data
        end
        seq = data isa Insertion ? data.seq : data
        encoded = _as_seq(S, seq)
        isempty(encoded) && throw(ArgumentError("insertion sequence must not be empty"))
        return Insertion(encoded)
    elseif kind == DeltaTypeSnp
        base = data isa Substitution ? data.base : data
        return Substitution(_as_symbol(S, base))
    elseif kind == DeltaTypeSV
        if data isa StructuralVariant{S}
            isempty(data.seq) && throw(ArgumentError("structural variant sequence must not be empty"))
            return data
        end
        seq = data isa StructuralVariant ? data.seq : data
        encoded = _as_seq(S, seq)
        isempty(encoded) && throw(ArgumentError("structural variant sequence must not be empty"))
        return StructuralVariant(encoded)
    elseif kind == DeltaTypeCNV
        if data isa CopyNumber{S}
            data.copies >= 1 || throw(ArgumentError("copy number must be positive, got $(data.copies)"))
            isempty(data.seq) && throw(ArgumentError("copy-number sequence must not be empty"))
            return data
        end
        if data isa CopyNumber
            seq, copies = data.seq, data.copies
        else
            seq, copies = data
        end
        copies >= 1 || throw(ArgumentError("copy number must be positive, got $copies"))
        encoded = _as_seq(S, seq)
        isempty(encoded) && throw(ArgumentError("copy-number sequence must not be empty"))
        return CopyNumber(encoded, Int(copies))
    else
        throw(ArgumentError("unsupported delta type $kind"))
    end
end

"""
    JournalEntry(kind, position, data, time)
    JournalEntry(kind, interval::ReferenceInterval, data, time)
    JournalEntry{S}(kind, position, data, time)

One edit against a reference of sequence type `S`.

`position` is a 1-based coordinate on the reference. Replay maps it through
every earlier edit, so removing an event leaves the other coordinates valid.
`position` must be at least 1.

`data` is normalized into a [`DeltaPayload`](@ref):

  - `DeltaTypeDel`: a positive length, or a `Deletion`
  - `DeltaTypeIns`: a sequence, or an `Insertion`
  - `DeltaTypeSnp`: a symbol, a `Char`, or a `Substitution`
  - `DeltaTypeSV`: a sequence, or a `StructuralVariant`
  - `DeltaTypeCNV`: `(sequence, copies)`, or a `CopyNumber`

Allele sequences are copied. A later mutation of the caller's sequence does
not change the stored edit.

The interval constructor reads the coordinate from a [`ReferenceInterval`](@ref).
A deletion uses a non-empty interval and ignores `data` as the length. A
substitution requires an interval of length 1. An insertion, copy-number
change, or structural variant requires an empty interval anchored before
`first`.

`time` orders edits and is the key in a [`DeltaMap`](@ref). A negative time
raises `DomainError`.

# Examples
```jldoctest
julia> entry = JournalEntry(DeltaTypeIns, 2, "TT", 0);

julia> String(entry.payload.seq)
"TT"

julia> isempty(locus(entry))
true
```
"""
struct JournalEntry{S <: LongSequence}
    delta_type :: DeltaType
    position :: Int64
    payload :: DeltaPayload{S}
    time :: Timestamp

    function JournalEntry{S}(kind :: DeltaType, position :: Integer, data, time :: Integer) where {S <: LongSequence}
        pos = Int64(position)
        pos >= 1 || throw(ArgumentError("reference position must be >= 1, got $pos"))
        stamped = time isa Timestamp ? time : Timestamp(time)
        # Copy allele sequences so a later edit of the caller's sequence
        # cannot rewrite history already stored in the journal.
        return new{S}(kind, pos, _copy_payload(make_payload(S, kind, data)), stamped)
    end
end

_copy_payload(payload :: Deletion) = payload
_copy_payload(payload :: Substitution) = payload
_copy_payload(payload :: Insertion) = Insertion(copy(payload.seq))
_copy_payload(payload :: StructuralVariant) = StructuralVariant(copy(payload.seq))
_copy_payload(payload :: CopyNumber) = CopyNumber(copy(payload.seq), payload.copies)

function JournalEntry(kind :: DeltaType, position :: Integer, data, time :: Integer)
    sequence_type = _sequence_type_of(data)
    return JournalEntry{sequence_type}(kind, position, data, time)
end

function JournalEntry{S}(kind :: DeltaType, interval :: ReferenceInterval, data, time :: Integer) where {S <: LongSequence}
    if kind == DeltaTypeDel
        isempty(interval) &&
            throw(ArgumentError("a deletion interval must contain at least one reference base"))
        return JournalEntry{S}(kind, first(interval), length(interval), time)
    elseif kind == DeltaTypeSnp
        length(interval) == 1 ||
            throw(ArgumentError("a substitution interval must cover exactly one reference base"))
        return JournalEntry{S}(kind, first(interval), data, time)
    else
        isempty(interval) ||
            throw(ArgumentError("an inserted allele uses an empty ReferenceInterval anchored before `first`"))
        return JournalEntry{S}(kind, first(interval), data, time)
    end
end

function JournalEntry(kind :: DeltaType, interval :: ReferenceInterval, data, time :: Integer)
    sequence_type = _sequence_type_of(data)
    return JournalEntry{sequence_type}(kind, interval, data, time)
end

"""
    locus(entry) -> ReferenceInterval

Reference span affected by `entry`, as a 1-based closed interval.

A deletion covers `position:position+len-1`. A substitution covers the
single base `position`. An insertion, copy-number change, or structural
variant is an empty interval anchored immediately before `position`.

# Examples
```jldoctest
julia> locus(JournalEntry(DeltaTypeDel, ReferenceInterval(2, 4), nothing, 0))
ReferenceInterval(2, 4)
```
"""
function locus(entry :: JournalEntry)
    position = entry.position
    payload = entry.payload
    if payload isa Deletion
        return ReferenceInterval(position, position + payload.len - 1)
    elseif payload isa Substitution
        return ReferenceInterval(position, position)
    else
        return ReferenceInterval(position, position - 1)
    end
end

# Sequences are mutable, so Julia does not synthesize content equality for
# the payload structs that hold them.
Base.:(==)(a :: Insertion, b :: Insertion) = a.seq == b.seq
Base.hash(a :: Insertion, h :: UInt) = hash(a.seq, hash(Insertion, h))

Base.:(==)(a :: StructuralVariant, b :: StructuralVariant) = a.seq == b.seq
Base.hash(a :: StructuralVariant, h :: UInt) = hash(a.seq, hash(StructuralVariant, h))

Base.:(==)(a :: CopyNumber, b :: CopyNumber) = a.copies == b.copies && a.seq == b.seq
Base.hash(a :: CopyNumber, h :: UInt) = hash(a.copies, hash(a.seq, hash(CopyNumber, h)))

Base.:(==)(a :: JournalEntry, b :: JournalEntry) =
    a.delta_type == b.delta_type &&
    a.position == b.position &&
    a.payload == b.payload &&
    a.time == b.time
Base.hash(a :: JournalEntry, h :: UInt) =
    hash(a.time, hash(a.payload, hash(a.position, hash(a.delta_type, h))))

Base.isless(a :: JournalEntry, b :: JournalEntry) = a.position < b.position

function Base.show(io :: IO, entry :: JournalEntry)
    print(io, "JournalEntry(", entry.delta_type,
          ", ref=", entry.position,
          ", ", entry.payload,
          ", t=", UInt64(entry.time), ")")
end


"""
    DeltaMap{S}()
    DeltaMap(n::Integer) -> Vector{DeltaMap}
    DeltaMap(n::Integer, ::Type{S}) -> Vector{DeltaMap{S}}

Time-ordered map from a [`Timestamp`](@ref) to a [`JournalEntry`](@ref) for
reference sequences of type `S`.

`DeltaMap(n, S)` builds `n` empty journals. `n` must be at least 1. The
one-argument form uses `LongDNA{4}`.
"""
const DeltaMap{S <: LongSequence} = SortedDict{Timestamp, JournalEntry{S}, Base.Order.ForwardOrdering}

function DeltaMap{S}() where {S <: LongSequence}
    return SortedDict{Timestamp, JournalEntry{S}}(Base.Order.ForwardOrdering())
end

DeltaMap(n :: Integer) = DeltaMap(n, LongDNA{4})
function DeltaMap(n :: Integer, ::Type{S}) where {S <: LongSequence}
    n >= 1 || throw(ArgumentError("a journal needs at least one sequence, got $n"))
    return [DeltaMap{S}() for _ in 1:n]
end
DeltaMap() = DeltaMap{LongDNA{4}}()


"""
    JournaledString(reference)
    JournaledString(reference, n::Integer)
    JournaledString(reference, maps, time=0, id="")

A reference sequence plus one journal of edits for each derived sequence.

`reference` may be any `LongSequence`, including `LongDNA{2}`, `LongDNA{4}`,
`LongRNA{2}`, and `LongRNA{4}`. The stored reference is a copy. Calling
`JournaledString(reference)` creates 10 derived sequences.
`JournaledString(reference, n)` creates `n` sequences, with `n >= 1`.
Passing `maps` keeps those per-sequence journals as they are, so a later
edit is visible on the vector the caller still holds.

`current_time` is the [`Timestamp`](@ref) assigned to the next edit.
`reference_id` is the FASTA or FASTQ identifier when the reference was
loaded from a file, and `""` otherwise.

Positions stored in the journals are 1-based reference coordinates.
[`sequence`](@ref) materializes a derived sequence.

# Examples
```jldoctest
julia> journal = JournaledString(dna"ACGT", 1);

julia> add_delta!(journal, 1, DeltaTypeIns, 2, "TT");

julia> String(sequence(journal, 1))
"ATTCGT"
```
"""
mutable struct JournaledString{S <: LongSequence}
    reference :: S
    deltaMap :: Vector{DeltaMap{S}}
    current_time :: Timestamp
    reference_id :: String
end

function JournaledString(
    reference :: S,
    maps :: AbstractVector,
    time :: Integer = 0,
    id :: AbstractString = "",
) where {S <: LongSequence}
    return JournaledString{S}(
        copy(reference),
        Vector{DeltaMap{S}}(maps),
        time isa Timestamp ? time : Timestamp(time),
        String(id),
    )
end

function JournaledString(reference :: S, n :: Integer) where {S <: LongSequence}
    return JournaledString(reference, DeltaMap(n, S))
end

JournaledString(reference :: LongSequence) = JournaledString(reference, 10)


"""
    JSTNode{S}

One node in the tree derived from a [`JSTree`](@ref) journal.

`parent` is `nothing` at a root. `delta` maps a sequence index to the
[`JournalEntry`](@ref) shared by that branch. `children` holds the nodes
that follow this one.
"""
mutable struct JSTNode{S <: LongSequence}
    parent :: Union{Nothing, JSTNode{S}}
    delta :: Dict{Int64, JournalEntry{S}}
    children :: Vector{JSTNode{S}}
end

"""
    JSTree(reference, n)

Journal of `n` derived sequences, with a tree kept for display.

`journal` is the source of truth: a map from a reference position to the
edit buckets at that position. `rootChildren` is rebuilt from `journal`
by `print_tree`. `is_equal` compares `root`, `current_time`, `length`,
and `journal`. `n` must be at least 1. The stored reference is a copy.

# Examples
```jldoctest
julia> tree = JSTree(dna"ACGT", 1);

julia> add_delta!(tree, 1, DeltaTypeSnp, 1, 'T');

julia> String(sequence(tree, 1))
"TCGT"
```
"""
mutable struct JSTree{S <: LongSequence}
    root :: S
    rootChildren :: Vector{JSTNode{S}}
    journal :: SortedDict{Int64, Vector{Dict{Int64, JournalEntry{S}}}}
    current_time :: Timestamp
    length :: Int
end

function JSTree(root :: S, n_sequences :: Integer) where {S <: LongSequence}
    n_sequences >= 1 || throw(ArgumentError("a tree needs at least one sequence, got $n_sequences"))
    return JSTree{S}(
        copy(root),
        JSTNode{S}[],
        SortedDict{Int64, Vector{Dict{Int64, JournalEntry{S}}}}(),
        Timestamp(0),
        Int(n_sequences),
    )
end

### struct.jl ends here.
