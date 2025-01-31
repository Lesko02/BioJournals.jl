### -*- Mode: Julia -*-

using BioJournals
using Test
using BioSequences
using DataStructures


@testset "BioJournals.jl" begin
    reference_seq = LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAG")
    deltaMap = [SortedDict{Int, JournalEntry}() for _ in 1:10]
    jst = JournaledString(reference_seq, deltaMap)
    add_delta!(jst, [1, 2], DeltaTypeIns, 8, "CGTA")
    add_delta!(jst, [4, 9], DeltaTypeSnp, 10, 'C')
    add_delta!(jst, [8], DeltaTypeIns, 24, "NNNNN")
    
    js1 = JournaledString(reference_seq, [SortedDict{Int, JournalEntry}() for _ in 1:10], 0)
    js2 = JournaledString(reference_seq, [SortedDict{Int, JournalEntry}() for _ in 1:10], 0)
    add_delta!(js1, [1, 2], DeltaTypeIns, 8, "CGTA")
    add_delta!(js2, [1, 2], DeltaTypeIns, 8, "CGTA")

    @test (js1==js2) == false

    @test (js1==js1) == true

    @test is_equal(js1, js2) == true
    
    @test Base.isequal(js1, js2) == false

    @test build_sequences(jst) ==
        "Sequence 1: AGATCGAATGCGCGAGCTAGCGACTCAG\n" *
        "Sequence 2: AGATCGAATGCGCGAGCTAGCGACTCAG\n" *
        "Sequence 3: AGATCGAGCGAGCTAGCGACTCAG\n" *
        "Sequence 4: AGATCGAGCCAGCTAGCGACTCAG\n" *
        "Sequence 5: AGATCGAGCGAGCTAGCGACTCAG\n" *
        "Sequence 6: AGATCGAGCGAGCTAGCGACTCAG\n" *
        "Sequence 7: AGATCGAGCGAGCTAGCGACTCAG\n" *
        "Sequence 8: AGATCGAGCGAGCTAGCGACTCANNNNNG\n" *
        "Sequence 9: AGATCGAGCCAGCTAGCGACTCAG\n" *
        "Sequence 10: AGATCGAGCGAGCTAGCGACTCAG\n"
    

end

### runtests.jl ends here.
