### -*- Mode: Julia -*-

### structs.jl

"""
Represents a time stamp tagging a given "change".

The time stamp is a primitive type, stored as an unsigned 64-bit integer.
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
Enumeration of delta types for sequence modifications.

# Values: 
   - `DeltaTypeDel`: Deletion 
   - `DeltaTypeIns`: Insertion 
   - `DeltaTypeSnp`: Single nucleotide polymorphism 
   - `DeltaTypeSV`: Structural variation 
   - `DeltaTypeCNV`: Copy number variation
"""
@enum DeltaType DeltaTypeDel DeltaTypeIns DeltaTypeSnp DeltaTypeSV DeltaTypeCNV

"""`DeltaTypeDel`: remove one or more reference bases."""
DeltaTypeDel

"""`DeltaTypeIns`: insert a sequence immediately before a reference base."""
DeltaTypeIns

"""`DeltaTypeSnp`: replace one reference symbol."""
DeltaTypeSnp

"""`DeltaTypeSV`: insert a structural-variant allele, or append it past the reference end."""
DeltaTypeSV

"""`DeltaTypeCNV`: insert tandem copies of a sequence before a reference base."""
DeltaTypeCNV

"""
    ReferenceInterval(first, last)

1-based closed interval on the reference. Both `first` and `last` are
included, the same convention as `GenomicFeatures.Interval`.

An empty interval (`last == first - 1`) is an insertion anchor: the new
bases sit immediately before reference base `first`. In VCF that event is
written with `POS == first - 1`, because VCF names the base before the allele.
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

Remove `len` reference bases starting at the entry's reference position.
"""
struct Deletion
    len :: Int
end

"""
    Insertion(seq)

Insert `seq` immediately before the reference base at the entry's position.
`seq` uses the same alphabet as the journal reference.
"""
struct Insertion{S <: LongSequence}
    seq :: S
end

"""
    Substitution(base)

Replace one reference symbol.
"""
struct Substitution{T}
    base :: T
end

"""
    StructuralVariant(seq)

Allele inserted by a structural-variant edit.

When the entry position names an existing reference base, `seq` is inserted
immediately before that base. A position at or beyond the last reference
base appends `seq` instead. The number of leading gap symbols is
`position - length(reference)`. A positive pad requires a 4-bit alphabet
(`LongDNA{4}` or `LongRNA{4}`). `LongDNA{2}` and `LongRNA{2}` have no gap.
"""
struct StructuralVariant{S <: LongSequence}
    seq :: S
end

"""
    CopyNumber(seq, copies)

Insert `copies` tandem copies of `seq` before the reference base at the
entry's position.
"""
struct CopyNumber{S <: LongSequence}
    seq :: S
    copies :: Int
end

"""Edit stored in a [`JournalEntry`](@ref) whose reference sequence type is `S`."""
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

One mutation against the reference sequence.

`position` is a 1-based coordinate on the **reference**, not on an
already edited copy. Replay maps it through every earlier edit, so
removing an event leaves the other reference coordinates valid.

`data` is normalized into a [`DeltaPayload`](@ref). `time` orders edits
and is the key used by [`DeltaMap`](@ref).
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

Reference span affected by `entry`, using a 1-based closed interval.

Deletions and substitutions cover the bases they change. Insertions, copy-number
changes, and interior structural variants are empty intervals anchored immediately
before the reference base named by `entry.position`.
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
Sorted mapping of time keys to journal entries for reference sequences of type `S`.
"""
const DeltaMap{S <: LongSequence} = SortedDict{Timestamp, JournalEntry{S}, Base.Order.ForwardOrdering}

function DeltaMap{S}() where {S <: LongSequence}
    return SortedDict{Timestamp, JournalEntry{S}}(Base.Order.ForwardOrdering())
end

"""
    DeltaMap(n) -> Vector
    DeltaMap(n, S) -> Vector

`n` empty journals whose reference sequence type is `S`. `S` defaults to
`LongDNA{4}`, the IUPAC DNA alphabet used by FASTA.
"""
DeltaMap(n :: Integer) = DeltaMap(n, LongDNA{4})
function DeltaMap(n :: Integer, ::Type{S}) where {S <: LongSequence}
    n >= 1 || throw(ArgumentError("a journal needs at least one sequence, got $n"))
    return [DeltaMap{S}() for _ in 1:n]
end
DeltaMap() = DeltaMap{LongDNA{4}}()


"""
A reference sequence plus one journal of edits per derived sequence.

`reference` may be any `LongSequence`, including `LongDNA{2}`, `LongDNA{4}`,
`LongRNA{2}`, and `LongRNA{4}`. The stored reference is a copy. Edits store
1-based reference coordinates. The per-sequence maps passed to the
constructor are kept as-is, so an edit made through the journal is visible
on that same vector. `reference_id` is the FASTA/FASTQ identifier when the
reference was loaded from a file, and `""` otherwise.
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
One node in the tree derived from a [`JSTree`](@ref) journal.
"""
mutable struct JSTNode{S <: LongSequence}
    parent :: Union{Nothing, JSTNode{S}}
    delta :: Dict{Int64, JournalEntry{S}}
    children :: Vector{JSTNode{S}}
end

"""
A journaled string tree. `journal` is the source of truth. `rootChildren`
is rebuilt from it by [`build_tree!`](@ref) for display.
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
