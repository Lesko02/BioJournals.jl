### -*- Mode: Julia -*-

### journal.jl

@testset "replay" begin
    reference = dna"ACGTACGT"

    @test apply_delta(reference, JournalEntry[]) == reference
    @test String(apply_delta(reference, JournalEntry(DeltaTypeSV, 2, dna"TT", 0))) == "ATTCGTACGT"
    @test String(apply_delta(reference, JournalEntry(DeltaTypeSV, 10, dna"TT", 0))) == "ACGTACGT--TT"
    @test String(apply_delta(reference, JournalEntry(DeltaTypeCNV, 2, (dna"GG", 1), 0))) == "AGGCGTACGT"

    @test_throws ArgumentError apply_delta(reference, JournalEntry(DeltaTypeSnp, 9, 'A', 0))
    @test_throws ArgumentError apply_delta(reference, JournalEntry(DeltaTypeDel, 7, 3, 0))
    @test_throws ArgumentError apply_delta(reference, [
        JournalEntry(DeltaTypeDel, 2, 2, 0),
        JournalEntry(DeltaTypeSnp, 3, 'A', 1),
    ])

    journal = JournaledString(reference)
    @test_throws BoundsError BioJournals.edits_at(journal, 0)
    tree = JSTree(reference, 2)
    @test_throws BoundsError BioJournals.edits_at(tree, 3)
    @test isempty(BioJournals.edits_at(tree, 1))
end

function capture_stdout(f)
    original = stdout
    reader, writer = redirect_stdout()
    try
        f()
        flush(stdout)
    finally
        redirect_stdout(original)
        close(writer)
    end
    return read(reader, String)
end

@testset "tree journal" begin
    reference = dna"ACGTACGTACGT"
    tree = JSTree(reference, 3)
    add_delta!(tree, [1], DeltaTypeIns, 4, "AA")
    add_delta!(tree, [2], DeltaTypeIns, 4, "TT")
    add_delta!(tree, [1], DeltaTypeSnp, 8, 'C')
    add_delta!(tree, [2], DeltaTypeSnp, 8, 'G')

    @test BioJournals.delta_count(tree) == 4
    @test length(BioJournals.edits_at(tree, 1)) == 2
    printed = capture_stdout() do
        print_tree(tree)
    end
    @test occursin("type=DeltaTypeIns", printed)
    @test occursin("type=DeltaTypeSnp", printed)
    @test count(==('\n'), printed) >= 3

    sequence_type = LongDNA{4}
    buckets = [
        Dict{Int64, JournalEntry{sequence_type}}(),
        Dict(Int64(1) => JournalEntry{sequence_type}(DeltaTypeSnp, 1, 'A', 0)),
    ]
    fresh = BioJournals._matching_bucket!(sequence_type, buckets, DeltaTypeIns, "TT")
    @test isempty(fresh)
    @test length(buckets) == 3

    tree.journal[99] = Dict{Int64, JournalEntry{LongDNA{4}}}[]
    BioJournals.build_tree!(tree)
    @test BioJournals.delta_count(tree) == 4

    alone = JSTree(reference, 1)
    add_delta!(alone, [1], DeltaTypeSnp, 1, 'A')
    add_delta!(alone, [1], DeltaTypeSnp, 4, 'T')
    remove_delta!(alone, 1)
    @test BioJournals.delta_count(alone) == 1
    remove_delta!(alone, 0)
    @test BioJournals.delta_count(alone) == 0
    @test isempty(alone.journal)
    @test_throws ArgumentError remove_delta!(alone, 4)

    history = JournaledString(dna"ACGT", DeltaMap(1))
    add_delta!(history, [1], DeltaTypeSnp, 2, 'A')
    @test occursin("Time 0", get_mutation_history(history.deltaMap[1]))
    @test occursin("at index 1", get_mutation_interval(history, 0, 0))
    @test get_mutation_interval(history, 5, 6) == ""
    @test String(only(get_sequences_at_time(history, 0))) == "AAGT"

    named = Dict(
        "root" => UnitRange{Int64}[],
        "leaf" => UnitRange{Int64}[1:3],
    )
    report = capture_stdout() do
        print_results(named)
    end
    @test occursin("No Match at root", report)
    @test occursin("Match at leaf", report)
    @test occursin("1:3", report)
end

### journal.jl ends here.
