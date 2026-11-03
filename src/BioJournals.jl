"""
    BioJournals

Journaled DNA and RNA sequences for mutation tracking.

Edits are stored as 1-based reference coordinates and replayed in time
order. [`sequence`](@ref) materializes a derived sequence as an ordinary
`LongSequence`.
"""
module BioJournals

using BioSequences
using DataStructures
using FASTX
using CodecZlib

import FASTX: sequence

include("structs.jl")
include("api.jl")
include("functions.jl")
include("apply.jl")
include("algorithms.jl")
include("io.jl")

export DeltaType,
    DeltaTypeDel,
    DeltaTypeIns,
    DeltaTypeSnp,
    DeltaTypeSV,
    DeltaTypeCNV,
    Deletion,
    Insertion,
    Substitution,
    StructuralVariant,
    CopyNumber,
    JournalEntry,
    JournaledString,
    JSTree,
    DeltaMap,
    Timestamp,
    ReferenceInterval,
    locus,
    sequence,
    sequences,
    nsequences,
    read_journal,
    edits_at,
    print_tree,
    print_sequences,
    build_sequences,
    print_deltas,
    add_delta!,
    remove_delta!,
    apply_delta,
    approximate_search,
    exact_search,
    get_mutation_history,
    get_mutation_interval,
    get_sequences_at_time,
    simulate_mutation!,
    remove_mutation!,
    is_equal,
    print_results,
    load_fasta,
    write_fasta,
    load_fastq,
    write_fastq

end
