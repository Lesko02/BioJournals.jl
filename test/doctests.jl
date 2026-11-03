using Documenter
using BioJournals
using BioSequences

DocMeta.setdocmeta!(
    BioJournals,
    :DocTestSetup,
    :(using BioJournals, BioSequences);
    recursive = true,
)

@testset "Docstrings" begin
    doctest(BioJournals; manual = false)
end
