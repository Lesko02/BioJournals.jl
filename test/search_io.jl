### -*- Mode: Julia -*-

### search_io.jl

@testset "search windows" begin
    reference = dna"ACGTACGT"
    inserted = JournaledString(reference, DeltaMap(1))
    add_delta!(inserted, [1], DeltaTypeIns, 2, "TT")
    needle = dna"TT"
    found = exact_search(inserted, needle)[1]
    materialized = apply_delta(reference, BioJournals.edits_at(inserted, 1))
    @test found == collect(findall(ExactSearchQuery(needle), materialized))
    @test found == [2:3]

    changed = JournaledString(dna"AGATCGAGCGAGCTAGCGACTCAG", DeltaMap(1))
    add_delta!(changed, [1], DeltaTypeSnp, 13, 'A')
    @test exact_search(changed, dna"GCTA")[1] == UnitRange{Int64}[]

    removed = JournaledString(dna"AGATCGAGCGAGCTAGCGACTCAG", DeltaMap(1))
    add_delta!(removed, [1], DeltaTypeDel, 20, 5)
    @test exact_search(removed, dna"TCAG")[1] == UnitRange{Int64}[]

    appended = JournaledString(dna"AGATCGAGCGAGCTAGCGACTCAG", DeltaMap(1))
    add_delta!(appended, [1], DeltaTypeSV, 26, "AA")
    tail = exact_search(appended, dna"TCAG")[1]
    rebuilt = apply_delta(appended.reference, BioJournals.edits_at(appended, 1))
    @test tail == collect(findall(ExactSearchQuery(dna"TCAG"), rebuilt))
    @test String(rebuilt) == "AGATCGAGCGAGCTAGCGACTCAG--AA"

    swallowed = JournaledString(dna"ACGTACGT", DeltaMap(1))
    add_delta!(swallowed, [1], DeltaTypeIns, 3, "TT")
    add_delta!(swallowed, [1], DeltaTypeDel, 2, 2)
    @test exact_search(swallowed, dna"AT")[1] == [1:2]

    cleared = JournaledString(dna"ACGT", DeltaMap(1))
    add_delta!(cleared, [1], DeltaTypeDel, 1, 4)
    @test exact_search(cleared, dna"AC")[1] == UnitRange{Int64}[]

    shifted = apply_delta(dna"ACGTACGT", [
        JournalEntry(DeltaTypeSV, 2, dna"TT", 0),
        JournalEntry(DeltaTypeSnp, 4, 'A', 1),
    ])
    @test String(shifted) == "ATTCGAACGT"
end

@testset "FASTQ validation" begin
    sequence = dna"ACGT"
    newline_id = Dict("bad\nid" => (sequence, "!!!!"))
    short_quality = Dict("read" => (sequence, "!!"))
    path = tempname()
    try
        @test_throws ErrorException write_fastq(newline_id, path)
        @test_throws ErrorException write_fastq(short_quality, path * ".gz")
    finally
        rm(path; force = true)
        rm(path * ".gz"; force = true)
    end
end

### search_io.jl ends here.
