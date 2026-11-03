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

Read FASTA records from `filename`.

DNA is returned as `LongDNA{4}`. A record whose sequence contains `U` and
no `T` is returned as `LongRNA{4}`. A name ending in `.gz` is decompressed.
A missing file raises an error. A repeated record identifier raises
`ArgumentError`.
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
    write_fasta(sequences, filename)
    write_fasta(journal, filename)

Write FASTA records to `filename`, wrapping lines at 70 columns.

`sequences` maps an identifier to a `LongSequence`. `journal` writes every
derived sequence in index order. The identifier is `reference_id_<index>`
when the reference was loaded from a file, and `sequence_<index>` otherwise.
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
    write_fastq(records, filename)

Write FASTQ records to `filename`.

`records` maps an identifier to `(sequence, quality)`. The quality string
must have the same number of code units as the sequence. An identifier
that contains a newline raises an error. A name ending in `.gz` is
compressed. The writer stores a quality header.
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
    load_fastq(filename) -> Dict{String, Tuple{LongSequence, String}}

Read FASTQ records from `filename`.

Each value is `(sequence, quality)`. DNA is `LongDNA{4}`. A record whose
sequence contains `U` and no `T` is `LongRNA{4}`. A name ending in `.gz`
is decompressed. A missing file raises an error. A repeated identifier
raises `ArgumentError`.
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

FASTQ quality is not kept. The format follows the file name: `.fasta`,
`.fa`, `.fastq`, or `.fq`, with an optional `.gz` suffix. `id` selects a
record by identifier. Otherwise record `index` is used, in file order and
1-based. `n` is the number of derived sequences, all unedited. A missing
file, an empty file, or an unknown `id` raises an error. `n` must be at
least 1.
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
