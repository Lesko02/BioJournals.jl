### -*- Mode: Julia -*-

### io.jl


"""
Loads sequences from a FASTA file into a dictionary.

# Args:
   - `filename`: Path to the FASTA file.

# Returns:
   - A dictionary mapping sequence IDs to `LongDNA{4}` sequences.

# Raises:
   - An error if the file does not exist.
"""
function load_fasta(filename)
    if !isfile(filename)
        error("File '$filename' does not exist.")
    end
    
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


"""
Writes sequences to a FASTA file.

Creates a new file if `filename` does not exist.

# Args:
   - `sequences`: A dictionary mapping sequence IDs to `LongDNA{4}` sequences.
   - `filename`: Path to the output FASTA file.
"""
function write_fasta(sequences :: Dict{String, LongDNA{4}}, filename :: String)
    open(filename, "w") do io
        writer = FASTA.Writer(io, 70) 
        for (id, seq) in sequences
            record = FASTA.Record(id, seq)
            write(writer, record)
        end
        close(writer)
    end
end


"""
Writes records to a FASTQ file.

Creates a new file if `filename` does not exist.

# Args:
   - `records`: A dictionary mapping sequence IDs to tuples of `LongDNA{4}`
                sequences and their quality scores.
   - `filename`: Path to the output FASTQ file.
"""
function write_fastq(records :: Dict{String, Tuple{LongDNA{4}, String}},
                     filename :: String)
    io = endswith(filename, ".gz") ? 
        GzipCompressorStream(open(filename, "w")) : 
        open(filename, "w")

    try
        writer = FASTQ.Writer(io, true)
        for (id, (seq, qual)) in records

            seq_str = String(seq) 
            qual_str = String(qual)

            if occursin(r"[\n\r]", id)
                error("ID contains newline characters: ", repr(id))
            end

            if length(seq_str) != length(qual_str)
                error("""
                        Length mismatch in record $(repr(id)):
                        - Sequence: $(length(seq_str)) bases
                        - Quality: $(length(qual_str)) chars
                      """)
            end

            record_data = Vector{UInt8}()
            append!(record_data, "@", id, "\n")
            append!(record_data, seq_str, "\n+\n")
            append!(record_data, qual_str, "\n")

            record = FASTQ.Record(record_data)
            write(writer, record)
        end
        close(writer)
    finally
        close(io)
    end
end


"""
Loads records from a FASTQ file into a dictionary.

# Args:
   - `filename`: Path to the FASTQ file.

# Returns:
   - A dictionary mapping sequence IDs to tuples of `LongDNA{4}` sequences 
     and their quality scores.

# Raises:
   - An error if the file does not exist.
"""
function load_fastq(filename :: String)
    if !isfile(filename)
        error("File '$filename' does not exist.")
    end
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


### io.jl ends here.
