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
Base.:(==)(ts::Timestamp, x::Integer) = x >= 0 && UInt64(ts) == UInt64(x)
Base.:(==)(x::Integer, ts::Timestamp) = ts == x
Base.:(<)(ts::Timestamp, x::Integer) = x >= 0 && UInt64(ts) < UInt64(x)
Base.:(<)(x::Integer, ts::Timestamp) = x < 0 || UInt64(x) < UInt64(ts)
Base.:<=(ts::Timestamp, x::Integer) = x >= 0 && UInt64(ts) <= UInt64(x)
Base.:<=(x::Integer, ts::Timestamp) = x < 0 || UInt64(x) <= UInt64(ts)

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
"""
struct Insertion
    seq :: LongDNA{4}
end

"""
    Substitution(base)

Replace the reference base at the entry's position.
"""
struct Substitution
    base :: DNA
end

"""
    StructuralVariant(seq)

Insert `seq` before a reference base, or append it when the position is at
or beyond the reference end. Positions past the end are padded with `-`.
"""
struct StructuralVariant
    seq :: LongDNA{4}
end

"""
    CopyNumber(seq, copies)

Insert `copies` tandem copies of `seq` before the reference base at the
entry's position.
"""
struct CopyNumber
    seq :: LongDNA{4}
    copies :: Int
end

"""Payload stored on a [`JournalEntry`](@ref). Positions are reference coordinates."""
const DeltaPayload = Union{Deletion, Insertion, Substitution, StructuralVariant, CopyNumber}

_as_dna4(seq :: LongDNA{4}) = copy(seq)
_as_dna4(seq) = LongDNA{4}(seq)
_as_base(base :: DNA) = base
_as_base(base) = convert(DNA, base)

"""
Translate a constructor argument into a [`DeltaPayload`](@ref).

Accepted inputs, by delta type:

  - deletion: a positive length, or a `Deletion`
  - insertion: a DNA sequence, or an `Insertion`
  - SNP: a `DNA` or `Char`, or a `Substitution`
  - structural variant: a DNA sequence, or a `StructuralVariant`
  - copy number: `(sequence, copies)`, or a `CopyNumber`
"""
function make_payload(kind :: DeltaType, data)
    if kind == DeltaTypeDel
        data isa Deletion && return data
        len = Int(data)
        len >= 1 || throw(ArgumentError("deletion length must be positive, got $len"))
        return Deletion(len)
    elseif kind == DeltaTypeIns
        data isa Insertion && return data
        seq = _as_dna4(data)
        isempty(seq) && throw(ArgumentError("insertion sequence must not be empty"))
        return Insertion(seq)
    elseif kind == DeltaTypeSnp
        data isa Substitution && return data
        return Substitution(_as_base(data))
    elseif kind == DeltaTypeSV
        data isa StructuralVariant && return data
        seq = _as_dna4(data)
        isempty(seq) && throw(ArgumentError("structural variant sequence must not be empty"))
        return StructuralVariant(seq)
    elseif kind == DeltaTypeCNV
        data isa CopyNumber && return data
        seq, copies = data
        copies >= 1 || throw(ArgumentError("copy number must be positive, got $copies"))
        return CopyNumber(_as_dna4(seq), Int(copies))
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
struct JournalEntry
    delta_type :: DeltaType
    position :: Int64
    payload :: DeltaPayload
    time :: Timestamp

    function JournalEntry(kind :: DeltaType, position :: Integer, data, time :: Integer)
        pos = Int64(position)
        pos >= 1 || throw(ArgumentError("reference position must be >= 1, got $pos"))
        stamped = time isa Timestamp ? time : Timestamp(time)
        return new(kind, pos, make_payload(kind, data), stamped)
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
Sorted mapping of time keys to journal entries.

# Keys: 
   - `Timestamp` time keys. 

# Values: 
   - Corresponding `JournalEntry` objects.
"""
const DeltaMap = SortedDict{Timestamp, JournalEntry, Base.Order.ForwardOrdering}

DeltaMap() = SortedDict{Timestamp, JournalEntry}(Base.Order.ForwardOrdering())

function DeltaMap(n :: Int)
    d = [SortedDict{Timestamp, JournalEntry}() for _ in 1:n]
    return d
end


"""
A reference DNA sequence plus one journal of edits per derived sequence.
Edits store 1-based reference coordinates.

# Fields: 
   - `reference`: Original LongDNA{4} sequence. 
   - `deltaMap`: Vector of DeltaMaps storing modifications. 
   - `current_time`: Current timestamp for modifications.

# Constructor: 
    JournaledString(reference, deltaMap)
   - Initializes with `current_time` set to 0.
"""
mutable struct JournaledString
    reference :: LongDNA{4}
    deltaMap  :: Vector{DeltaMap}
    current_time :: Timestamp
end


function JournaledString(reference :: LongDNA{4}, deltaMap :: Vector{DeltaMap})
    JournaledString(reference, deltaMap, 0)
end


function JournaledString(reference :: LongDNA{4})
    JournaledString(reference, DeltaMap(10))
end



"""
A node in a Journaled String Tree (JST).

# Fields: 
   - `parent`: Reference to the parent node (or nothing if root). 
   - `delta`: Modifications applied at this node. 
   - `children`: Children of the node.
"""
mutable struct JSTNode
    parent   :: Union{Nothing, JSTNode}
    delta    :: Dict{Int64, JournalEntry}
    children :: Vector{JSTNode}
end

"""
A Journaled String Tree (JST).

# Fields: 
   - `root`: The root sequence (LongDNA{4}). 
   - `children`: Dictionary mapping node names to JSTNode objects.

# Constructor: 
   - JSTree(root_sequence) 
   - Initializes with a `root` node named "root".
"""
mutable struct JSTree
    root         :: LongDNA{4}
    rootChildren :: Vector{JSTNode}
    journal      :: SortedDict{Int64, Vector{Dict{Int64, JournalEntry}}}
    current_time :: Timestamp
    length       :: Int
end

function JSTree(root::LongDNA{4}, default_length::Int)
    return JSTree(
        root,
        JSTNode[],
        SortedDict{Int64, Vector{Dict{Int64,JournalEntry}}}(),
        0,
        default_length
    )
end

### struct.jl ends here.
