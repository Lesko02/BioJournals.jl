### -*- Mode: Julia -*-

### runtests.jl

using Test
using BioSequences
using DataStructures
using BioJournals

@testset "BioJournals.jl" begin
    reference_seq = LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAG")
    deltaMap = DeltaMap(10)
    jst = JournaledString(reference_seq, deltaMap)
    add_delta!(jst, [1, 2], DeltaTypeIns, 8, "CGTA")
    add_delta!(jst, [4, 9], DeltaTypeSnp, 10, 'C')
    add_delta!(jst, [8], DeltaTypeIns, 24, "NNNNN")

    ## Test of apply_delta
    
    entry1 = JournalEntry(DeltaTypeDel, 1, 24, 1)
    entry2 = JournalEntry(DeltaTypeSV, 24, LongDNA{4}("NN") , 2)
    entry4 = JournalEntry(DeltaTypeCNV, 1, (LongDNA{4}("ATCG"), 2), 3)
    @test apply_delta(reference_seq, entry1) == LongDNA{4}("")
    @test apply_delta(reference_seq, entry2) == 
        LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAGNN")
    @test apply_delta(reference_seq, entry4) == 
        LongDNA{4}("ATCGATCGAGATCGAGCGAGCTAGCGACTCAG") 
    @test reference_seq == apply_delta(reference_seq, deltaMap[10])
    @test reference_seq != apply_delta(reference_seq, deltaMap[1])
    @test reference_seq != apply_delta(reference_seq, deltaMap[2])

    ## end Test of apply_delta

    ## Test of comparisons JSS

    js1 = JournaledString(reference_seq,
    DeltaMap(10), 0)
    js2 = JournaledString(reference_seq, DeltaMap(10), 0)
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


    ## end Test of comparisons JSS

    ## Test of comparisons JST
    
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

    ## end Test of comparisons JST

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

    ## Done like this to avoid withespace problems

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


    ## Test of flatten
    
    @test flatten(tree1, "child1") == LongDNA{4}("AGATCGACGTAGCGAGCTAGCGACTCAG")

    ## end Test of flatten

    ## Test of removal
    
    @test length(tree1.children) == 6
    remove_node!(tree1, "child1")
    @test length(tree1.children) == 1

    ## end Test of removal


    ## Test of search
    
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
    
    @test approximate_search(js1, needle2, 5) == test_return3
    @test approximate_search(js1, needle2) == test_return3

    ## Tree section
    
    test_return4 = Dict{String, Vector{UnitRange{Int64}}}(
        name => UnitRange{Int64}[] for name in keys(tree2.children))
    
    voidvector = UnitRange{Int64}[]
    test_return5 = Dict( 
        "child1" => [16:19],
        "root"   => voidvector,
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

    ## end Test of search


    ## Testing for Errors
    
    @test_throws ErrorException add_node(tree2, "child18", deltaMap[2] , "child2")
    @test_throws ErrorException remove_node!(tree1, "child18")
    @test_throws ErrorException remove_node!(tree1, "root")
    @test_throws ErrorException approximate_search(js1, needle, 115) == test_return

    ## end Testing for Errors

    ## "Testing" the prints to update the code coverage
    js3 = JournaledString(LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAG"),
                          DeltaMap(10),
                          0)

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
        print_tree(tree2)
        print_deltas(js3)
        get_mutation_history(js3.deltaMap[1])
        get_mutation_interval(js3, 1, 5)
        print_results(test_return)
        print_results(test_return2)
        print_results(test_return5)
        print_results(test_return6)
    end

    ## End "Testing" the prints to update the code coverage

    ## Test of trim_node
    
    @test length(tree2.children) == 6
    @test_throws ErrorException trim_node(tree2, "child18")
    @test_throws ErrorException trim_node(tree2, "root")
    trim_node(tree2, "child4")
    @test length(tree2.children) == 5

    ## End Test of trim_node

    ## Test of IO functions
    
    ioTestDict = Dict("Try1" => reference_seq)
    mktemp() do path, io
        write_fasta(ioTestDict, path)
        @test_throws ErrorException load_fasta("InvalidPATH.!/")
        close(io)
        ioTestDict2 = load_fasta(path)
        @test ioTestDict == ioTestDict2
    end

    ioTestDict3 = Dict("Try2" => (reference_seq, "!A5@DECC>GID<7?BEGC?=:7!"))
    mktemp() do path, io
        write_fastq(ioTestDict3, path)
        @test_throws ErrorException load_fastq("InvalidPATH.!/")
        close(io)
        ioTestDict4 = load_fastq(path)
        @test [x[1] for (id, x) in ioTestDict3] == [x[1] for (id, x) in ioTestDict4]
    end
    
    ## End Test of IO functions
    
    ## Testing new apis
    
    @test_throws ErrorException add_delta!(tree2, "root", entry3)
    @test_throws ErrorException add_delta!(tree2, "nonesiste", entry3)
    @test_throws ErrorException add_delta!(tree2, "root", DeltaTypeSnp, 21, 'C')
    @test_throws ErrorException add_delta!(tree2, "nonesiste", DeltaTypeSnp, 21, 'C')
    add_delta!(tree2, "child1", DeltaTypeSnp, 21, 'C')
    add_delta!(tree2, "child1", entry3)
    @test length(tree2.children["child1"].deltaMap) == 3
    println(tree2.children["child1"])
    @test_throws ErrorException remove_delta!(tree2, "child1", 3)
    remove_delta!(tree2, "child1", 2)
    @test length(tree2.children["child1"].deltaMap) == 2
    add_node(tree2, "root", DeltaMap(), "child6")

    ## end testing new apis
end

test_entries = [
    # Deletion: Remove 2 bases starting at position 3
    JournalEntry(DeltaTypeDel, 3, 2, 1),
    
    # Insertion: Add "AC" at position 2
    JournalEntry(DeltaTypeIns, 2, "AC", 2),
    
    # SNP: Change base at position 4 to 'T'
    JournalEntry(DeltaTypeSnp, 4, DNA_T, 3),
    
    # Structural Variation: Add "TT" at position 5
    JournalEntry(DeltaTypeSV, 5, "TT", 4),
    
    # Copy Number Variation: Repeat "GG" 3 times at position 6
    JournalEntry(DeltaTypeCNV, 6, (LongDNA{4}("GG"), 3), 5),

    JournalEntry(DeltaTypeCNV, 6, (LongDNA{4}("GG"), 3), 5),
    
    JournalEntry(DeltaTypeCNV, 6, (LongDNA{4}("GG"), 3), 5)
]
tree14 = JSTree2(LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAGCCCCCC"))

for entry in test_entries
    add_delta2!(tree14, entry)
end

print_tree2(tree14)
### runtests.jl ends here.
