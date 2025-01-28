# Definition of DeltaTypes
@enum DeltaType DeltaTypeDel DeltaTypeIns DeltaTypeSnp DeltaTypeSV DeltaTypeCNV

# Journal Entry
struct JournalEntry
    delta_type::DeltaType    # Type of delta
    position::Int            # Index of delta
    data::Any                # Additional Data
    time::Int                # Timestamp
end

# Definition of JournaledString Structure
mutable struct JournaledString
    reference::LongDNA{4}
    deltaMap::Vector{SortedDict{Int, JournalEntry, Base.Order.ForwardOrdering}}
    current_time::Int
end

# Constructor for JournaledString
function JournaledString(reference::LongDNA{4},
    deltaMap::Vector{SortedDict{Int, JournalEntry, Base.Order.ForwardOrdering}})
    JournaledString(reference, deltaMap, 0)  # Default value for current_time
end

#JST
struct JSTNode
    parent::Union{Nothing, JSTNode}
    deltaMap::Union{Nothing, SortedDict{Int, JournalEntry, Base.Order.ForwardOrdering}}
    name::String
end

struct JSTree
    root::LongDNA{4};
    children::Dict{String, JSTNode}
 end

 function JSTree(root_sequence::LongDNA{4})
    # Create the root node with "Nothing" as the parent and an empty deltaMap
    root_node = JSTNode(nothing, nothing, "root")
    return JSTree(root_sequence, Dict("root" => root_node))
end

function add_node(tree::JSTree, parent_name::String, 
    deltas::SortedDict{Int, JournalEntry, Base.Order.ForwardOrdering}, 
    node_name::String)

    if !haskey(tree.children, parent_name)
        error("Parent node '$parent_name' does not exist.")
    else
        # If parent is not root, ensure the parent exists in the tree
        parent_node = tree.children[parent_name]
        new_node = JSTNode(parent_node, deltas, node_name)
    end

# Add the new node to the tree
tree.children[node_name] = new_node
end