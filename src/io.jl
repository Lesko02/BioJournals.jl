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

function write_fastq(records::Dict{String, Tuple{LongDNA{4}, String}},
    filename::String)
    io = endswith(filename, ".gz") ? 
         GzipCompressorStream(open(filename, "w")) : 
         open(filename, "w")

    try
        writer = FASTQ.Writer(io, true) 
        for (id, (seq, qual)) in records
            record = FASTQ.Record(id, nothing, string(seq), collect(qual))
            write(writer, record)
        end
        close(writer)
    finally
        close(io)
    end
end



function load_fastq(filename::String)
    records = Dict{String, Tuple{LongDNA{4}, String}}()
    io = endswith(filename, ".gz") ? 
    GzipDecompressorStream(open(filename, "r")) : open(filename, "r")
    try
        reader = FASTQ.Reader(io)
        for record in reader
            id = String(FASTX.identifier(record))
            seq = LongDNA{4}(FASTX.sequence(record))
            qual = String(FASTX.quality(record))
            records[id] = (seq, qual)
        end
        close(reader)
    finally
        close(io)
    end
    return records
end
