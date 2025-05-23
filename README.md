# BioJournals

[![Julia](https://img.shields.io/badge/Julia-1.11.1-blue)](https://julialang.org/)
[![Build Status](https://github.com/Lesko02/BioJournals.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Lesko02/BioJournals.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Lesko02/BioJournals.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Lesko02/BioJournals.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Project Status: WIP](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

Simulators tracking mutations in a cell population, e.g., to study the
evolution of a solid tumor, must manipulate many DNA/RNA sequences at
a time.  It is obviously unfeasible to duplicate a full sequence every
time a mutation happens.  Any such simulator must therefore rely
on sophisticated data structures just tracking the *delta*'s from one
sequence to the next, while keeping them in a *tree-like* (or more
complex) data structure.


## Goal

The goal of the the `BioJournals` project is to provide a **Julia**
library that implements a form of *journaled sequence*, alongside a
proper API for its manipulation.  An important *key* identifying such
sequences will be *time*.

Similar data structures that have been proposed to solve similar
problems (cf., [SeqAn](https://www.seqan.de/data-structures/)), but
they target different problems that the one at hand.  Plus, we would
like a **Julia** implementation to complement our other efforts in the
(solid tumors) simulation arena.


## Installation

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


## Usage

Create and manipulate journaled sequences:

```julia  
using BioJournals

reference_seq = LongDNA{4}("AGATCGAGCGAGCTAGCGACTCAG")            # Creating a reference sequence
deltaMap = DeltaMap(10)                                           # Creating a DeltaMap object
jst = JournaledString(reference_seq, deltaMap)                    # Inizialization

add_delta!(jst, [1], DeltaTypeIns, 8, "CGTA")                     # Adding deltas such as Insertions
add_delta!(jst, [2], DeltaTypeSnp, 10, 'C')                       # and Permutations

print_sequences(jst)                                              # Printing the modified sequences
```

## Documentation

Get interactive help using the package directly from the REPL, use:

```julia
?JournaledString    # To view the docstring about the JournaledString type
```

produces this in the REPL:

```console
search: JournaledString JournalEntry unsafe_string

  A DNA sequence with associated modifications.

  Fields:
  ≡≡≡≡≡≡≡

    •  reference: Original LongDNA{4} sequence.

    •  deltaMap: Vector of DeltaMaps storing modifications.

    •  current_time: Current timestamp for modifications.

  Constructor:
  ≡≡≡≡≡≡≡≡≡≡≡≡

  JournaledString(reference, deltaMap)

    •  Initializes with current_time set to 0.
```
