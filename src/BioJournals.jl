### -*- Mode: Julia -*-

### BioJournals.jl

module BioJournals

using BioSequences
using DataStructures
using FASTX
using CodecZlib

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
    JSTNode,

    DeltaMap,
    Timestamp,
    
    insert!,
    delete_at!,
    structure_variation!,
    copy_number_variation!,
    
    print_tree,
    print_sequences,
    build_sequences,
    print_deltas,
    
    add_delta!,
    remove_delta!,
    apply_delta,
    
    approximate_findall,
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
    write_fastq,

    build_tree!,
    delta_count
    
end

### BioJournls.jl ends here.
