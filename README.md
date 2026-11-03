# BioJournals

[![Julia](https://img.shields.io/badge/Julia-1.12%20%7C%201.13-blue)](https://julialang.org/)
[![Build Status](https://github.com/Lesko02/BioJournals.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Lesko02/BioJournals.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Lesko02/BioJournals.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Lesko02/BioJournals.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Project Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

Simulators tracking mutations in a cell population, e.g., to study the
evolution of a solid tumor, must manipulate many DNA/RNA sequences at
a time.  It is obviously unfeasible to duplicate a full sequence every
time a mutation happens.  Any such simulator must therefore rely
on sophisticated data structures just tracking the *delta*'s from one
sequence to the next, while keeping them in a *tree-like* (or more
complex) data structure.


## Goal

The goal of the `BioJournals` project is to provide a **Julia**
library that implements a form of *journaled sequence*, alongside a
proper API for its manipulation.  An important *key* identifying such
sequences will be *time*.

Similar data structures that have been proposed to solve similar
problems (cf., [SeqAn](https://www.seqan.de/data-structures/)), but
they target different problems that the one at hand.  Plus, we would
like a **Julia** implementation to complement our other efforts in the
(solid tumors) simulation arena.


## Installation

BioJournals requires Julia 1.12 or newer.

The current release is [v1.0.0](https://github.com/Lesko02/BioJournals.jl/releases/tag/v1.0.0).

```julia
using Pkg
Pkg.add(url="https://github.com/Lesko02/BioJournals.jl", rev="v1.0.0")
```

From the Pkg REPL (`]`):

```julia
add https://github.com/Lesko02/BioJournals.jl#v1.0.0
```

Omit the revision to follow the default branch.

The package builds on [BioSequences](https://github.com/BioJulia/BioSequences.jl) for alphabets and search, and on [FASTX](https://github.com/BioJulia/FASTX.jl) for FASTA and FASTQ.


## Quick start

A journal keeps one reference and a list of derived sequences. Edits are stored against the reference. `sequence` replays them into an ordinary `LongSequence`.

```julia
using BioJournals
using BioSequences

reference = LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAG")
journal = JournaledString(reference, 2)

add_delta!(journal, 1, DeltaTypeIns, 8, "CGTA")   # before reference base 8
add_delta!(journal, 2, DeltaTypeSnp, 10, 'C')     # reference base 10, G -> C

String(sequence(journal, 1))
# "AGATCGACGTAGCGAGCTAGCGACTCAG"

String(sequence(journal, 2))
# "AGATCGAGCCAGCTAGCGACTCAG"
```

The same two edits on one sequence still use reference coordinates. The insertion does not move the SNP:

```julia
stacked = JournaledString(reference, 1)
add_delta!(stacked, 1, DeltaTypeIns, 8, "CGTA")
add_delta!(stacked, 1, DeltaTypeSnp, 10, 'C')
String(sequence(stacked, 1))
# "AGATCGACGTAGCCAGCTAGCGACTCAG"
```

`sequence` is `FASTX.sequence`, so `using BioJournals, FASTX` is one function. The returned sequence can be passed to `findall`, `translate`, and the FASTX writers.

Every call to `add_delta!` stamps the journal clock and then advances it by one. The first edit is time `0`. `remove_delta!(journal, index, time)` drops that edit and leaves the clock where it is, so the next edit receives a new time. The other reference coordinates stay valid.


## Reference coordinates

Positions are 1-based indexes on the **reference**, in the same closed-interval convention as `GenomicFeatures`.

| Edit | Where the allele goes |
|:--|:--|
| Insertion, copy number, interior structural variant | immediately before reference base `position` |
| Substitution | replaces reference base `position` |
| Deletion of length `len` | removes `position:position+len-1` |
| Structural variant with `position >= length(reference)` | appended after the reference |

A structural variant at or past the end writes `position - length(reference)` gap symbols before the allele. `position == length(reference)` therefore appends with no gap. `position == length(reference) + 2` writes two gaps, then the allele.

`ReferenceInterval(first, last)` names that span. Both ends are included. `ReferenceInterval(2, 4)` covers bases 2, 3, and 4. An empty interval (`last == first - 1`) is an insertion anchor before base `first`:

```julia
locus(JournalEntry(DeltaTypeIns, 8, "CGTA", 0))
# ReferenceInterval(8, before)

locus(JournalEntry(DeltaTypeDel, ReferenceInterval(2, 4), nothing, 0))
# ReferenceInterval(2, 4)
```

In VCF, that insertion is written with `POS == first - 1`, because VCF names the base before the allele. `JournalEntry(kind, interval, data, time)` accepts the interval directly. A deletion requires a non-empty interval. A substitution requires length 1. An inserted allele requires the empty anchor.


## Alphabets

`JournaledString` and `JSTree` accept `LongDNA{2}`, `LongDNA{4}`, `LongRNA{2}`, and `LongRNA{4}`. The stored reference and every allele are copies, so later changes to the caller's sequence do not rewrite the journal.

| Alphabet | Symbols | Gap |
|:--|:--|:--|
| `LongDNA{4}`, `LongRNA{4}` | IUPAC, including `N` | yes |
| `LongDNA{2}`, `LongRNA{2}` | `ACGT` / `ACGU` | no |

Recording a substitution outside the alphabet, or a structural variant that needs a gap on a 2-bit alphabet, raises `ArgumentError` at `add_delta!`. FASTA and FASTQ readers build a 4-bit sequence: DNA, or RNA when the record contains `U` and no `T`.

```julia
rna = JournaledString(LongRNA{4}("ACGUACGU"), 1)
add_delta!(rna, 1, DeltaTypeSnp, 4, 'C')
String(sequence(rna, 1))
# "ACGCACGU"
```


## Edits

`add_delta!(journal, index, kind, position, data)` records one edit. `index` may be a single integer or a vector of 1-based sequence indexes. A vector stores the same allele once per index, at successive timestamps.

| `kind` | `data` |
|:--|:--|
| `DeltaTypeDel` | positive length, or a `Deletion` |
| `DeltaTypeIns` | sequence string or `LongSequence`, or an `Insertion` |
| `DeltaTypeSnp` | `Char`, a symbol of the alphabet, or a `Substitution` |
| `DeltaTypeSV` | sequence, or a `StructuralVariant` |
| `DeltaTypeCNV` | `(sequence, copies)`, or a `CopyNumber` |

```julia
String(apply_delta(dna"ACGT", JournalEntry(DeltaTypeSV, 6, dna"TT", 0)))
# "ACGT--TT"

String(apply_delta(dna"ACGT", JournalEntry(DeltaTypeCNV, 2, (dna"GG", 3), 0)))
# "AGGGGGGCGT"
```

`edits_at(journal, i)` returns the edits of one sequence, oldest first. `get_sequences_at_time(journal, t)` replays only edits at or before `t`. `is_equal` compares the reference, the clock, and the journals. For a string it also compares `reference_id`.


## Search

`exact_search` and `approximate_search` scan every derived sequence without rebuilding the ones that have no edits. Hits are 1-based closed intervals in the derived sequence.

```julia
edited = JournaledString(dna"ACGTACGT", 1)
add_delta!(edited, 1, DeltaTypeIns, 2, "TT")
exact_search(edited, dna"TT")[1]
# 1-element Vector{UnitRange{Int64}}:
#  2:3
```

`approximate_search(journal, needle, tol)` allows `tol` percent mismatches, from 1 through 99. The default is 5. The error budget is `ceil(length(needle) * tol / 100)`. After a hit, the next search starts `budget + 1` bases later, matching BioSequences. `print_results` writes the dictionary to `stdout`.


## Files

```julia
cohort = read_journal("tumour.fasta"; id="sample1", n=4)
add_delta!(cohort, 1, DeltaTypeSnp, 120, 'A')
write_fasta(cohort, "derived.fasta")   # sample1_1, sample1_2, ...
```

`read_journal` also accepts FASTQ (`.fastq`, `.fq`) and gzip (`.gz`). Quality scores are not stored. `index` selects a record in file order when `id` is omitted.

`load_fasta` and `load_fastq` return dictionaries. `write_fasta` and `write_fastq` accept those dictionaries. A repeated record identifier raises `ArgumentError`. FASTQ quality must have the same number of code units as the sequence, and an identifier cannot contain a newline.


## Trees

`JSTree(reference, n)` stores the same journal as `JournaledString`. `print_tree` rebuilds a branch view: edits that share a position, kind, and allele share a branch. `sequence`, `add_delta!`, `remove_delta!`, and the search functions accept a tree.


## API

The names exported by `using BioJournals` are the supported interface. In the REPL, `?JournaledString` prints the manual entry.

| Topic | Names |
|:--|:--|
| Clock and intervals | `Timestamp`, `ReferenceInterval`, `locus` |
| Edits | `DeltaType`, `DeltaTypeDel`, `DeltaTypeIns`, `DeltaTypeSnp`, `DeltaTypeSV`, `DeltaTypeCNV`, `Deletion`, `Insertion`, `Substitution`, `StructuralVariant`, `CopyNumber`, `JournalEntry` |
| Journals | `JournaledString`, `JSTree`, `DeltaMap`, `add_delta!`, `remove_delta!`, `simulate_mutation!`, `remove_mutation!`, `edits_at`, `is_equal` |
| Sequences | `sequence`, `sequences`, `nsequences`, `apply_delta`, `get_sequences_at_time` |
| Reports | `print_sequences`, `build_sequences`, `print_deltas`, `print_tree`, `get_mutation_history`, `get_mutation_interval`, `print_results` |
| Search | `exact_search`, `approximate_search` |
| Files | `read_journal`, `load_fasta`, `write_fasta`, `load_fastq`, `write_fastq` |


## Testing

```julia
using Pkg
Pkg.test("BioJournals")
```

The suite includes [Aqua](https://github.com/JuliaTesting/Aqua.jl) checks for ambiguities, piracy, compat bounds, and undocumented public names. Examples in the docstrings are checked with [Documenter](https://github.com/JuliaDocs/Documenter.jl).


## Related packages

- [BioSequences.jl](https://github.com/BioJulia/BioSequences.jl) defines the alphabets, `LongSequence`, and the search queries used here.
- [FASTX.jl](https://github.com/BioJulia/FASTX.jl) reads and writes FASTA and FASTQ. `sequence` extends its generic function.
- [GenomicFeatures.jl](https://github.com/BioJulia/GenomicFeatures.jl) uses the same 1-based closed intervals. BioJournals does not depend on it.


## License

MIT. See [LICENSE](LICENSE).
