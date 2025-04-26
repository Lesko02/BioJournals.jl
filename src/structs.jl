### -*- Mode: Julia -*-

### structs.jl


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
    time :: Int                
end


"""
Sorted mapping of time keys to journal entries.

# Keys: 
   - `Integer` time keys. 

# Values: 
   - Corresponding `JournalEntry` objects.
"""
const DeltaMap = SortedDict{Int, JournalEntry, Base.Order.ForwardOrdering}

DeltaMap() = SortedDict{Int, JournalEntry}(Base.Order.ForwardOrdering())

function DeltaMap(n :: Int)
    d = [SortedDict{Int, JournalEntry}() for _ in 1:n]
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
    current_time :: Int
end


function JournaledString(reference :: LongDNA{4}, deltaMap :: Vector{DeltaMap})
    JournaledString(reference, deltaMap, 0)
end


"""
A node in a Journaled String Tree (JST).

# Fields: 
   - `parent`: Reference to the parent node (or nothing if root). 
   - `deltaMap`: Modifications applied at this node (or nothing). 
   - `name`: Name of the node.
"""
mutable struct JSTNode
    parent   :: Union{Nothing, JSTNode}
    deltaMap :: Union{Nothing, DeltaMap}
    name     :: String
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
    root::LongDNA{4}
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
function add_node(tree :: JSTree,
                  parent_name :: String,
                  deltas :: DeltaMap, 
                  node_name :: String)

    if !haskey(tree.children, parent_name)
        error("Parent node '$parent_name' does not exist.")
    else
        parent_node = tree.children[parent_name]
        new_node = JSTNode(parent_node, deltas, node_name)
    end
    tree.children[node_name] = new_node
end


"""
Removes a node and all its descendants from a JSTree.

# Arguments
- `tree::JSTree`: The journaled string tree.
- `node_name::String`: The name of the node to remove.

# Errors
- Throws an error if the node does not exist or if attempting to remove the root.
"""
function remove_node!(tree :: JSTree, node_name :: String)
    if node_name == "root"
        error("Cannot remove the root node.")
    end
    if !haskey(tree.children, node_name)
        error("Node '$node_name' does not exist.")
    end
    
    for (child_name, child_node) in collect(tree.children)
        if child_node.parent !== nothing && child_node.parent.name == node_name
            remove_node!(tree, child_name)
        end
    end
    delete!(tree.children, node_name)
end


function get_parent(jst :: JSTree, node_name :: String)
    node = jst.children[node_name]
    return node.parent
end


function trim_node(jst :: JSTree, node_name :: String)
    if node_name == "root"
        error("Cannot trim the root node.")
    end
    if !haskey(jst.children, node_name)
        error("Node '$node_name' does not exist.")
    end
    
    for (child_name, child_node) in collect(jst.children)
        if child_node.parent !== nothing && child_node.parent.name == node_name
            child_node.parent = get_parent(jst, node_name)
        end
    end
    delete!(jst.children, node_name)
end

### struct.jl ends here.
