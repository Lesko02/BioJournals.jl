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
    reader = FASTA.Reader(open(filename))
    try
        for record in reader
            id = String(FASTX.identifier(record))
            seq = LongDNA{4}(FASTX.sequence(String, record))
            seqs[id] = seq
        end
    finally
        close(reader)
    end
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
    writer = FASTA.Writer(open(filename, "w"); width = 70)
    try
        for (id, seq) in sequences
            write(writer, FASTA.Record(id, String(seq)))
        end
    finally
        close(writer)
    end
    return nothing
end


"""
Writes records to a FASTQ file.

Creates a new file if `filename` does not exist. A filename ending in `.gz`
is written through a gzip stream.

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
    writer = FASTQ.Writer(io; quality_header = true)
    try
        for (id, (seq, qual)) in records
            if occursin(r"[\n\r]", id)
                error("ID contains newline characters: ", repr(id))
            end

            seq_str = String(seq)
            qual_str = String(qual)
            if ncodeunits(seq_str) != ncodeunits(qual_str)
                error("""
                        Length mismatch in record $(repr(id)):
                        - Sequence: $(ncodeunits(seq_str)) bases
                        - Quality: $(ncodeunits(qual_str)) chars
                      """)
            end
            write(writer, FASTQ.Record(id, seq_str, qual_str))
        end
    finally
        close(writer)
    end
    return nothing
end


"""
Loads records from a FASTQ file into a dictionary.

A filename ending in `.gz` is read through a gzip stream.

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
        GzipDecompressorStream(open(filename, "r")) :
        open(filename, "r")
    reader = FASTQ.Reader(io)
    try
        for record in reader
            id = String(FASTX.identifier(record))
            seq = LongDNA{4}(FASTX.sequence(String, record))
            qual = String(FASTX.quality(record))
            records[id] = (seq, qual)
        end
    finally
        close(reader)
    end
    return records
end


### io.jl ends here.
