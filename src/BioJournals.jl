module BioJournals

using BioSequences
using DataStructures
using FASTX

include("structs.jl")
include("api.jl")
include("functions.jl")
include("algorithms.jl")

export DeltaType, DeltaTypeDel, DeltaTypeIns, DeltaTypeSnp, DeltaTypeSV, 
DeltaTypeCNV, JournalEntry, JournaledString, JSTree, JSTNode, DeltaMap,
insert!, delete_at!, structure_variation!, copy_number_variation!,
flatten, print_tree, print_sequences, build_sequences, print_deltas,
add_delta!, remove_delta!, apply_delta, add_node, remove_node!,
approximate_findall, approximate_search, exact_search,
get_mutation_history, get_mutation_interval, get_sequences_at_time,
simulate_mutation!, remove_mutation!, is_equal, print_results, trim_node
end
