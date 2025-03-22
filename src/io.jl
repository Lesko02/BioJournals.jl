function load_fasta(filename)
    seqs = Dict{String, LongDNA{4}}()
    reader = open(FASTA.Reader, filename)
    for record in reader
        id = String(FASTX.identifier(record))
        seq = LongDNA{4}(FASTX.sequence(record))
        seqs[id] = seq
    end
    close(reader)
    return seqs
end
#=
function write_fasta()
    # YET TO BE CODED
end
=#