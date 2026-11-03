### -*- Mode: Julia -*-

function _guess_sequence(text :: AbstractString)
    letters = uppercase(text)
    if occursin('U', letters) && !occursin('T', letters)
        return LongRNA{4}(text)
    else
        return LongDNA{4}(text)
    end
end

function _open_read(filename :: AbstractString)
    if endswith(lowercase(filename), ".gz")
        return GzipDecompressorStream(open(filename))
    else
        return open(filename)
    end
end

function _is_fastq(filename :: AbstractString)
    name = replace(lowercase(filename), r"\.gz$" => "")
    return endswith(name, ".fastq") || endswith(name, ".fq")
end

function _push_unique!(records, seen, id, value, filename)
    if id in seen
        throw(ArgumentError("duplicate record id $(repr(id)) in '$filename'"))
    end
    push!(seen, id)
    push!(records, (id, value))
    return records
end

"""
    load_fasta(filename) -> Dict{String, LongSequence}

Load FASTA records. DNA is read as `LongDNA{4}`. A record whose sequence
contains `U` and no `T` is read as `LongRNA{4}`. `.gz` files are decompressed.
"""
function load_fasta(filename)
    if !isfile(filename)
        error("File '$filename' does not exist.")
    end

    seqs = Dict{String, LongSequence}()
    seen = Set{String}()
    reader = FASTA.Reader(_open_read(filename))
    try
        for record in reader
            id = String(FASTX.identifier(record))
            id in seen && throw(ArgumentError("duplicate record id $(repr(id)) in '$filename'"))
            push!(seen, id)
            seqs[id] = _guess_sequence(FASTX.sequence(String, record))
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
function write_fasta(sequences :: AbstractDict{<:AbstractString, <:LongSequence}, filename :: AbstractString)
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
function write_fastq(records :: AbstractDict{<:AbstractString, <:Tuple{<:LongSequence, <:AbstractString}},
                     filename :: AbstractString)
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
    records = Dict{String, Tuple{LongSequence, String}}()
    seen = Set{String}()
    reader = FASTQ.Reader(_open_read(filename))
    try
        for record in reader
            id = String(FASTX.identifier(record))
            id in seen && throw(ArgumentError("duplicate record id $(repr(id)) in '$filename'"))
            push!(seen, id)
            seq = _guess_sequence(FASTX.sequence(String, record))
            qual = String(FASTX.quality(record))
            records[id] = (seq, qual)
        end
    finally
        close(reader)
    end
    return records
end

function _records_in_order(filename :: AbstractString)
    records = Tuple{String, LongSequence}[]
    seen = Set{String}()
    reader = _is_fastq(filename) ? FASTQ.Reader(_open_read(filename)) : FASTA.Reader(_open_read(filename))
    try
        for record in reader
            id = String(FASTX.identifier(record))
            _push_unique!(records, seen, id, _guess_sequence(FASTX.sequence(String, record)), filename)
        end
    finally
        close(reader)
    end
    return records
end

"""
    read_journal(filename; index=1, id=nothing, n=1) -> JournaledString

Build a journal whose reference is one record from a FASTA or FASTQ file.
FASTQ quality is not stored. The format is taken from the file name
(`.fasta`, `.fa`, `.fastq`, `.fq`, optionally compressed with `.gz`).

`id` selects a record by identifier. Otherwise record `index` is used in
file order (1-based). `n` is the number of derived sequences, all starting
unedited.
"""
function read_journal(filename :: AbstractString;
                      index :: Integer = 1,
                      id :: Union{Nothing, AbstractString} = nothing,
                      n :: Integer = 1)
    isfile(filename) || error("File '$filename' does not exist.")
    records = _records_in_order(filename)
    isempty(records) && throw(ArgumentError("file '$filename' contains no sequence records"))
    chosen = if id !== nothing
        found = findfirst(record -> record[1] == id, records)
        found === nothing && throw(ArgumentError("record id $(repr(id)) was not found in '$filename'"))
        records[found]
    else
        1 <= index <= length(records) || throw(BoundsError(records, index))
        records[index]
    end
    chosen_id, reference = chosen
    return JournaledString(reference, DeltaMap(n, typeof(reference)), 0, chosen_id)
end

"""
    write_fasta(journal, filename)

Write every derived sequence of `journal` as a FASTA record.
Identifiers are `reference_id` plus the sequence index when the reference
came from a file, and `sequence_<index>` otherwise.
"""
function write_fasta(journal :: JournaledString, filename :: AbstractString)
    prefix = isempty(journal.reference_id) ? "sequence" : journal.reference_id
    writer = FASTA.Writer(open(filename, "w"); width = 70)
    try
        for i in eachindex(journal.deltaMap)
            identifier = string(prefix, "_", i)
            write(writer, FASTA.Record(identifier, String(sequence(journal, i))))
        end
    finally
        close(writer)
    end
    return nothing
end

### io.jl ends here.
