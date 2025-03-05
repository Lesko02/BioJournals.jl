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
    delta_type::DeltaType    
    position::Int64          
    data::Any                
    time::Int                
end

"""
Sorted mapping of time keys to journal entries.

# Keys: 
   - `Integer` time keys. 

# Values: 
   - Corresponding `JournalEntry` objects.
"""
const DeltaMap = SortedDict{Int, JournalEntry, Base.Order.ForwardOrdering}

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
    reference::LongDNA{4}
    deltaMap::Vector{DeltaMap}
    current_time::Int
end

function JournaledString(reference::LongDNA{4}, deltaMap::Vector{DeltaMap})
    JournaledString(reference, deltaMap, 0)
end

"""
A node in a Journaled String Tree (JST).

# Fields: 
   - `parent`: Reference to the parent node (or nothing if root). 
   - `deltaMap`: Modifications applied at this node (or nothing). 
   - `name`: Name of the node.
"""
struct JSTNode
    parent::Union{Nothing, JSTNode}
    deltaMap::Union{Nothing, DeltaMap}
    name::String
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
struct JSTree
    root::LongDNA{4};
    children::Dict{String, JSTNode}
 end

 function JSTree(root_sequence::LongDNA{4})
    root_node = JSTNode(nothing, nothing, "root")
    return JSTree(root_sequence, Dict("root" => root_node))
end

"""
Adds a new node to a JSTree.

# Args: 
   - `tree`: The JSTree structure. 
   - `parent_name`: Name of the parent node. 
   - `deltas`: DeltaMap of modifications for the new node. 
   - `node_name`: Name of the new node.

# Returns: 
   - None (modifies the JSTree).

# Errors: 
   - Throws an error if the `parent` node does not exist.
"""
function add_node(tree::JSTree, parent_name::String, deltas::DeltaMap, 
    node_name::String)

        if !haskey(tree.children, parent_name)
            error("Parent node '$parent_name' does not exist.")
        else
            parent_node = tree.children[parent_name]
            new_node = JSTNode(parent_node, deltas, node_name)
        end
    tree.children[node_name] = new_node
end
