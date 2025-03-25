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

function write_fasta(sequences::Dict{String, LongDNA{4}}, filename::String)
    open(filename, "w") do io
        writer = FASTA.Writer(io, 70) 
        for (id, seq) in sequences
            record = FASTA.Record(id, seq)
            write(writer, record)
        end
        close(writer)
    end
end