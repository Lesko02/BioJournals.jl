### -*- Mode: Julia -*-

### sequences.jl

@testset "sequence primitives" begin
    @test String(BioJournals.insert!(LongDNA{4}("ACGT"), 2, dna"CG")) == "ACGCGT"
    @test String(delete_at!(LongDNA{4}("ACGT"), 2:1)) == "ACGT"
    @test String(delete_at!(LongDNA{4}("ACGTACGT"), 3:5)) == "ACCGT"

    @test String(structure_variation!(LongDNA{4}("ACGT"), 4, dna"TT")) == "ACGTTT"
    @test String(structure_variation!(LongDNA{4}("ACGT"), 2, dna"TT")) == "ATTCGT"
    @test String(structure_variation!(LongDNA{4}("ACGT"), 6, dna"TT")) == "ACGT--TT"

    @test String(copy_number_variation!(LongDNA{4}("ACGT"), 2, (dna"GG", 3))) == "AGGGGGGCGT"
end

@testset "payloads" begin
    inserted = Insertion(dna"AC")
    @test BioJournals.make_payload(DeltaTypeIns, inserted) == inserted
    @test BioJournals.make_payload(DeltaTypeIns, "AC") == inserted
    @test BioJournals.make_payload(DeltaTypeIns, dna"AC") == inserted
    @test hash(inserted) == hash(Insertion(dna"AC"))

    deleted = Deletion(2)
    @test BioJournals.make_payload(DeltaTypeDel, deleted) == deleted

    substituted = Substitution(DNA_A)
    @test BioJournals.make_payload(DeltaTypeSnp, substituted) == substituted
    @test BioJournals.make_payload(DeltaTypeSnp, DNA_C) == Substitution(DNA_C)

    variant = StructuralVariant(dna"TT")
    @test variant == StructuralVariant(dna"TT")
    @test BioJournals.make_payload(DeltaTypeSV, variant) === variant
    @test BioJournals.make_payload(DeltaTypeSV, dna"TT") == variant
    @test hash(variant) == hash(StructuralVariant(dna"TT"))

    copies = CopyNumber(dna"AT", 2)
    @test copies == CopyNumber(dna"AT", 2)
    @test BioJournals.make_payload(DeltaTypeCNV, copies) === copies
    @test BioJournals.make_payload(DeltaTypeCNV, (dna"AT", 2)) == copies
    @test hash(copies) == hash(CopyNumber(dna"AT", 2))

    left = JournalEntry(DeltaTypeIns, 2, inserted, Timestamp(1))
    right = JournalEntry(DeltaTypeIns, 2, "AC", 1)
    @test left == right
    @test hash(left) == hash(right)
    @test isless(left, JournalEntry(DeltaTypeIns, 3, "A", 0))
    @test occursin("ref=2", sprint(show, left))

    @test_throws ArgumentError BioJournals.make_payload(DeltaTypeIns, "")
    @test_throws ArgumentError BioJournals.make_payload(DeltaTypeDel, 0)
    @test_throws ArgumentError BioJournals.make_payload(DeltaTypeSV, dna"")
    @test_throws ArgumentError BioJournals.make_payload(DeltaTypeCNV, (dna"A", 0))
end

### sequences.jl ends here.
