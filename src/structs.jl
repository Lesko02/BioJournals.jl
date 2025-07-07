### -*- Mode: Julia -*-

### structs.jl

"""
Represents a time stamp tagging a given "change".

The time stamp is a `primitive` type, `Unisigned 64`.
"""
primitive type Timestamp <: Unsigned 64 end

Timestamp(x :: UInt8)  = reinterpret(Timestamp, x)
Timestamp(x :: UInt16) = reinterpret(Timestamp, x)
Timestamp(x :: UInt32) = reinterpret(Timestamp, x)
Timestamp(x :: UInt64) = reinterpret(Timestamp, x)
Timestamp(x::Unsigned) = reinterpret(Timestamp, x)
Timestamp(x::Int) = reinterpret(Timestamp, x)

UInt64(ts::Timestamp) = reinterpret(UInt64, ts)
Base.convert(::Type{Int64}, ts::Timestamp) = reinterpret(Int64, ts)
addto(ts::Timestamp, x::Unsigned) = Timestamp(UInt64(ts) + x)
addto(ts::Timestamp, x::Timestamp) = Timestamp(UInt64(ts) + UInt64(x))

Base.:+(ts::Timestamp, n::Unsigned) = Timestamp(UInt64(ts) + n)
Base.:+(ts::Timestamp, n::UInt64) = Timestamp(UInt64(ts) + n)

Base.:(==)(a::Timestamp, b::Timestamp) = UInt64(a) == UInt64(b)
Base.hash(ts::Timestamp, h::UInt64) = hash(UInt64(ts), h)

Base.:(<)(a::Timestamp, b::Timestamp) = UInt64(a) < UInt64(b)
Base.:(==)(ts::Timestamp, x::Integer) = UInt64(ts) == convert(UInt64, x)
Base.:(==)(x::Integer, ts::Timestamp) = ts == x
Base.:(<)(ts::Timestamp, x::Integer) = UInt64(ts) < convert(UInt64, x)
Base.:(<)(x::Integer, ts::Timestamp) = convert(UInt64, x) < UInt64(ts)
Base.:<=(ts::Timestamp, x::Integer) = UInt64(ts) <= convert(UInt64, x)
Base.:<=(x::Integer, ts::Timestamp) = convert(UInt64, x) <= UInt64(ts)

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
Represents a single modification in a journaled sequence.

# Fields: 
   - `delta_type`: Type of modification (DeltaType). 
   - `position`: Position of the modification. 
   - `data`: Modification data (can be any type). 
   - `time`: Time key for ordering modifications.
"""
struct JournalEntry
    delta_type :: DeltaType
    position :: Int64
    data :: Any
    time :: Timestamp
end

Base.isless(a::JournalEntry, b::JournalEntry) = a.position < b.position


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
A DNA sequence with associated modifications.

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
    journal = SortedDict{Int64, Dict{Int64, JournalEntry}}()
    return JSTree(
        root,
        JSTNode[],
        SortedDict{Int64, Vector{Dict{Int64,JournalEntry}}}(),
        0,
        default_length
    )
end

### struct.jl ends here.
