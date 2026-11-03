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

BioJournals requires Julia 1.12 or newer. Julia 1.10 cannot load the current
FASTX and BioSequences dependency set.

To install this package directly from GitHub, open the Julia REPL and 
run:  

```julia  
using Pkg  
Pkg.add(url="https://github.com/Lesko02/BioJournals.jl")  
```

Alternatively, use the Pkg REPL mode by pressing `]` and entering:

```julia
add https://github.com/Lesko02/BioJournals.jl
```


## Coordinates

Every mutation records a 1-based position on the **reference** sequence.
Replay maps that position through the edits that happened earlier, so an
insertion does not invalidate a later SNP, and deleting an edit does not
require rebasing the events that remain.

`ReferenceInterval(first, last)` is a 1-based closed interval: both ends
are included, the same convention as `GenomicFeatures`. An empty interval
(`last == first - 1`) anchors an insertion immediately before reference
base `first`. A VCF record for that insertion uses `POS == first - 1`.

An insertion, copy-number change, or structural variant at `position`
places its allele immediately before that reference base. A structural
variant whose position is at or beyond the last reference base is appended
instead. The number of leading gap symbols is `position - length(reference)`.
`LongDNA{2}` and `LongRNA{2}` cannot encode a gap, so a padded append on
those alphabets raises `ArgumentError`. Use `LongDNA{4}` or `LongRNA{4}`
when the variant runs past the reference.

## Usage

```julia
using BioJournals
using BioSequences

reference = LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAG")
journal = JournaledString(reference, 2)

add_delta!(journal, 1, DeltaTypeIns, 8, "CGTA")
add_delta!(journal, 2, DeltaTypeSnp, 10, 'C')

sequence(journal, 1)          # LongDNA{4}, usable by findall and translate
locus(edits_at(journal, 1)[1])

# A FASTA or FASTQ record can be the reference directly.
# cohort = read_journal("cohort.fasta"; id = "sample1", n = 4)
```

`sequence` is `FASTX.sequence`, so `using BioJournals, FASTX` does not make
the name ambiguous. `JournaledString` accepts `LongDNA{2}`, `LongDNA{4}`,
`LongRNA{2}`, and `LongRNA{4}`. FASTA and FASTQ readers build a 4-bit
sequence: DNA, unless the record contains `U` and no `T`, in which case RNA.

```julia
?JournaledString
?ReferenceInterval
?sequence
```
