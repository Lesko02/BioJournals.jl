using Aqua
using BioJournals

@testset "Aqua" begin
    Aqua.test_all(BioJournals; undocumented_names = true)
end
