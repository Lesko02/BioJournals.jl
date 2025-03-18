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

    jst2 = deepcopy(jst)
    jst2.reference = LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAGCCCCCC")

    @test (js1==js2) == false
    @test (js1==js1) == true
    @test is_equal(js1, js2) == true
    @test Base.isequal(js1, js2) == false
    @test is_equal(jst, js1) == false
    @test is_equal(jst, jst2) == false


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

    treecpy = JSTree(LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAGCCCCCC"))
    add_node(tree1, "root", deltaMap[1], "child1")
    add_node(tree1, "child1", deltaMap[2] , "child2")
    add_node(tree1, "child1", deltaMap[3], "child3")
    add_node(tree1, "child1", deltaMap[4], "child4")
    add_node(tree1, "child4", deltaMap[5], "child5")

    @test is_equal(tree1, tree2) == true
    @test Base.isequal(tree1, tree2) == false
    @test is_equal(tree1, treecpy) == false

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

    test_return_time = Vector{LongDNA{4}}()
    test_return_time = [
        dna"AGATCGACGTAGCGAGCTAGCGACTCAG",
        dna"AGATCGACGTAGCGAGCTAGCGACTCAG",
        dna"AGATCGAGCGAGCTAGCGACTCAG",
        dna"AGATCGAGCGAGCTAGCGACTCAGNNNNN",
        dna"AGATCGAGCGAGCTAGCGACTCAG",
        dna"AGATCGAGCGAGCTAGCGACTCAG",
        dna"AGATCGAGCGAGCTAGCGACTCAG",
        dna"AGATCGAGCGAGCTAGCGACTCAG",
        dna"AGATCGAGCGAGCTAGCGACTCAG",
        dna"AGATCGAGCGAGCTAGCGACCCAG"
    ]


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
    @test approximate_search(js1, needle, 15) == test_return
    
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

    test_return6 = Dict(
        "child1" => [16:18, 13:15, 20:23, 9:11],
        "root" => [16:18, 8:11, 13:15, 20:23, 9:11],
        "child2" => [9:11, 13:15, 24:27, 20:22, 16:18, 20:23],
        "child3" => [20:23, 8:11, 13:15, 24:27, 20:22, 9:11, 16:18],
        "child5" => [20:23, 8:11, 13:15, 24:27, 20:22, 9:11, 16:18],
        "child4" => [13:15, 24:27, 20:22, 20:23, 16:18]
    )

    @test exact_search(tree2, needle) == test_return4

    @test exact_search(tree2, needle2) == test_return5

    @test approximate_search(tree2, needle2) == test_return6

# end Test of search

# Testing for Errors
@test_throws ErrorException add_node(tree2, "child18", deltaMap[2] , "child2")
@test_throws ErrorException remove_node!(tree1, "child18")
@test_throws ErrorException remove_node!(tree1, "root")
@test_throws ErrorException approximate_search(js1, needle, 115) == test_return

# end Testing for Errors

# "Testing" the prints to update the code coverage
js3 = JournaledString(LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAG"),
[SortedDict{Int, JournalEntry}() for _ in 1:10], 0)

add_delta!(js3, [1, 2], DeltaTypeIns, 8, "CGTA")
add_delta!(js3, [10], DeltaTypeSnp, 21, 'C')
add_delta!(js3, [4], DeltaTypeSV, 24, dna"NNNNN")
add_delta!(js3, [1, 3], DeltaTypeCNV, 1, (LongDNA{4}("ATCG"), 2))
add_delta!(js3, [5, 6, 7], DeltaTypeIns, 5, "GTC")
add_delta!(js3, [8, 9], DeltaTypeDel, 10, 2)
add_delta!(js3, [4, 2], DeltaTypeSnp, 18, 'T')
add_delta!(js3, [1, 3], DeltaTypeSnp, 18, 'G')
add_delta!(js3, [6, 8], DeltaTypeSV, 20, dna"CCTG")
add_delta!(js3, [5, 7, 10], DeltaTypeDel, 3, 2)
add_delta!(js3, [4], DeltaTypeSV, 35, dna"ACTGACTG")
entry3 = JournalEntry(DeltaTypeIns, 3, "AGT", 10)
add_delta!(js3, [3, 4], entry3)
simulate_mutation!(js3, 5, entry3)
remove_mutation!(js3, 5, 23)
remove_delta!(js3, 3, 5)
@test get_sequences_at_time(js3, 3) == test_return_time
@test_throws ErrorException remove_delta!(js3, 1, 15)

redirect_stdout(devnull) do
    print_sequences(js3)
    print_sequences(tree2)
    print_tree(tree1)
    print_deltas(js3)
    get_mutation_history(js3.deltaMap[1])
    get_mutation_interval(js3.deltaMap[1], 1, 2)
    print_results(test_return2)
    print_results(test_return6)
end
# end "Testing" the prints to update the code coverage
end

### runtests.jl ends here.
