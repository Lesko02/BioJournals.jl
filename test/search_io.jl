### -*- Mode: Julia -*-

### search_io.jl

using FASTX

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

@testset "alphabets and intervals" begin
    dna2 = LongDNA{2}("ACGTACGT")
    narrow = JournaledString(dna2, 1)
    add_delta!(narrow, [1], DeltaTypeIns, 2, "TT")
    @test String(sequence(narrow, 1)) == "ATTCGTACGT"
    @test_throws ArgumentError add_delta!(narrow, [1], DeltaTypeSV, 10, "AA")

    rna = LongRNA{4}("ACGUACGU")
    rna_journal = JournaledString(rna, 1)
    add_delta!(rna_journal, [1], DeltaTypeSnp, 4, 'C')
    @test String(sequence(rna_journal, 1)) == "ACGCACGU"
    @test eltype(sequence(rna_journal, 1)) === RNA

    deletion = JournalEntry(DeltaTypeDel, ReferenceInterval(2, 4), nothing, 0)
    @test deletion.position == 2
    @test deletion.payload == Deletion(3)
    @test locus(deletion) == ReferenceInterval(2, 4)
    @test length(locus(deletion)) == 3
    @test 3 in locus(deletion)
    insertion = JournalEntry(DeltaTypeIns, ReferenceInterval(2, 1), "TT", 0)
    @test isempty(locus(insertion))
    @test first(locus(insertion)) == 2
    @test_throws ArgumentError JournalEntry(DeltaTypeDel, ReferenceInterval(2, 1), nothing, 0)
    @test_throws ArgumentError ReferenceInterval(0, 2)

    path = tempname() * ".fasta"
    try
        write_fasta(Dict("ref" => dna"ACGTACGT", "other" => dna"TTTTTTTT"), path)
        loaded = read_journal(path; id = "ref", n = 2)
        @test loaded.reference_id == "ref"
        @test nsequences(loaded) == 2
        @test String(loaded.reference) == "ACGTACGT"
        @test String(sequence(loaded, 1)) == "ACGTACGT"
        add_delta!(loaded, [1], DeltaTypeIns, 1, "GG")
        write_fasta(loaded, path)
        roundtrip = load_fasta(path)
        @test String(roundtrip["ref_1"]) == "GGACGTACGT"
        @test_throws ArgumentError read_journal(path; id = "missing")
    finally
        rm(path; force = true)
    end

    fq = tempname() * ".fastq"
    try
        write_fastq(Dict("read1" => (dna"ACGT", "IIII")), fq)
        from_fastq = read_journal(fq)
        @test from_fastq.reference_id == "read1"
        @test String(from_fastq.reference) == "ACGT"
    finally
        rm(fq; force = true)
    end
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

@testset "journal contracts" begin
    @test sequence === FASTX.sequence

    allele = dna"TT"
    owned = JournaledString(dna"ACGT", 1)
    add_delta!(owned, 1, DeltaTypeIns, 2, allele)
    push!(allele, DNA_A)
    @test String(sequence(owned, 1)) == "ATTCGT"

    ordered = JournaledString(dna"AC", 2)
    add_delta!(ordered, Int32[2], DeltaTypeIns, 1, "T")
    path = tempname() * ".fasta"
    try
        write_fasta(ordered, path)
        records = BioJournals._records_in_order(path)
        @test records[1][1] == "sequence_1"
        @test String(records[1][2]) == "AC"
        @test records[2][1] == "sequence_2"
        @test String(records[2][2]) == "TAC"
    finally
        rm(path; force = true)
    end

    narrow = JournaledString(LongDNA{2}("ACGT"), 1)
    @test_throws ArgumentError add_delta!(narrow, 1, DeltaTypeSnp, 2, 'N')
    rna2 = JournaledString(LongRNA{2}("ACGU"), 1)
    @test_throws ArgumentError add_delta!(rna2, 1, DeltaTypeSV, 6, "AA")
    rna4 = JournaledString(LongRNA{4}("ACGU"), 1)
    add_delta!(rna4, 1, DeltaTypeSV, 6, "AA")
    @test String(sequence(rna4, 1)) == "ACGU--AA"

    duplicate = tempname() * ".fasta"
    try
        open(duplicate, "w") do io
            write(io, ">a\nACGT\n>a\nTTTT\n")
        end
        @test_throws ArgumentError load_fasta(duplicate)
        @test_throws ArgumentError read_journal(duplicate)
    finally
        rm(duplicate; force = true)
    end
end

### search_io.jl ends here.
