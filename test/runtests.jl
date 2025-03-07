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

# Test of apply_delta
    entry1 = JournalEntry(DeltaTypeDel, 1, 24, 1)
    entry2 = JournalEntry(DeltaTypeSV, 24, LongDNA{4}("NN") , 2)
    @test apply_delta(reference_seq, entry1) == LongDNA{4}("")
    @test apply_delta(reference_seq, entry2) == LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAGNN")
    @test reference_seq == apply_delta(reference_seq, deltaMap[10])
    @test reference_seq != apply_delta(reference_seq, deltaMap[1])
    @test reference_seq != apply_delta(reference_seq, deltaMap[2])
# end Test of apply_delta

# Test of comparisons JSS
    js1 = JournaledString(reference_seq,
    [SortedDict{Int, JournalEntry}() for _ in 1:10], 0)
    js2 = JournaledString(reference_seq,
    [SortedDict{Int, JournalEntry}() for _ in 1:10], 0)
    add_delta!(js1, [1, 2], DeltaTypeIns, 8, "CGTA")
    add_delta!(js2, [1, 2], DeltaTypeIns, 8, "CGTA")

    @test (js1==js2) == false
    @test (js1==js1) == true
    @test is_equal(js1, js2) == true
    @test Base.isequal(js1, js2) == false

# end Test of comparisons JSS

# Test of comparisons JST
    tree1 = JSTree(LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAG"))
    add_node(tree1, "root", deltaMap[1], "child1")
    add_node(tree1, "child1", deltaMap[2] , "child2")
    add_node(tree1, "child1", deltaMap[3], "child3")
    add_node(tree1, "child1", deltaMap[4], "child4")
    add_node(tree1, "child4", deltaMap[5], "child5")

    tree2 = JSTree(LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAG"))
    add_node(tree2, "root", deltaMap[1], "child1")
    add_node(tree2, "child1", deltaMap[2] , "child2")
    add_node(tree2, "child1", deltaMap[3], "child3")
    add_node(tree2, "child1", deltaMap[4], "child4")
    add_node(tree2, "child4", deltaMap[5], "child5")

    @test is_equal(tree1, tree2) == true
    @test Base.isequal(tree1, tree2) == false
#end Test of comparisons JST

    test_return_build = "Sequence 1: AGATCGACGTAGCGAGCTAGCGACTCAG
Sequence 2: AGATCGACGTAGCGAGCTAGCGACTCAG
Sequence 3: AGATCGAGCGAGCTAGCGACTCAG
Sequence 4: AGATCGAGCCAGCTAGCGACTCAG
Sequence 5: AGATCGAGCGAGCTAGCGACTCAG
Sequence 6: AGATCGAGCGAGCTAGCGACTCAG
Sequence 7: AGATCGAGCGAGCTAGCGACTCAG
Sequence 8: AGATCGAGCGAGCTAGCGACTCANNNNNG
Sequence 9: AGATCGAGCCAGCTAGCGACTCAG
Sequence 10: AGATCGAGCGAGCTAGCGACTCAG"

    @test strip(build_sequences(jst)) == strip(test_return_build)
    #done like this to avoid withespace problems

# Test of flatten
    @test flatten(tree1, "child1") == LongDNA{4}("AGATCGACGTAGCGAGCTAGCGACTCAG")
# end Test of flatten

# Test of removal
    @test length(tree1.children) == 6
    remove_node!(tree1, "child1")
    @test length(tree1.children) == 1
# end Test of removal


# Test of search
    test_return = Dict(i => UnitRange{Int64}[] for i in 1:length(js1.deltaMap))
    test_return2= Dict(
        1 => [16:19],
        2 => [16:19],
        3 => [12:15],
        4 => [12:15],
        5 => [12:15],
        6 => [12:15],
        7 => [12:15],
        8 => [12:15],
        9 => [12:15],
        10 => [12:15]
    )

    test_return3= Dict(
        1 => [16:18, 13:15, 20:23, 9:11],
        2 => [16:18, 13:15, 20:23, 9:11],
        3 => [8:11, 13:15],
        4 => [8:11, 13:15],
        5 => [8:11, 13:15],
        6 => [8:11, 13:15],
        7 => [8:11, 13:15],
        8 => [8:11, 13:15],
        9 => [8:11, 13:15],
        10 => [8:11, 13:15]
    )
    
    needle = LongDNA{4}("AAAAAAAAAA")
    @test exact_search(js1, needle) == test_return

    needle2 = LongDNA{4}("GCTA")
    @test exact_search(js1, needle2) == test_return2
    

    @test approximate_search(js1, needle) == test_return
    
    @test approximate_search(js1, needle2) == test_return3

    #tree section
    test_return4 = Dict{String, Vector{UnitRange{Int64}}}(
        name => UnitRange{Int64}[] for name in keys(tree2.children))

    test_return5= Dict(
        "child1" => [16:19],
        "root" => [],
        "child2" => [20:23],
        "child3" => [16:19],
        "child5" => [16:19],
        "child4" => [16:19]
    )

    @test exact_search(tree2, needle) == test_return4

    @test exact_search(tree2, needle2) == test_return5

# end Test of search
end

### runtests.jl ends here.
