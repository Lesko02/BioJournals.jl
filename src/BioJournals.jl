module BioJournals

using BioSequences
using DataStructures
using FASTX

include("structs.jl")
include("functions.jl")

export insert!, delete_at!, structure_variation!, copy_number_variation!,
 JournalEntry, JournaledString, add_delta!, remove_delta!, apply_delta,
 print_sequences, build_sequences, print_deltas, get_mutation_history,
 get_mutation_interval, get_sequences_at_time, simulate_mutation!,
 remove_mutation!, is_equal
end
