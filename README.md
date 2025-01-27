# BioJournals

[![Build Status](https://github.com/Lesko02/BioJournals.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Lesko02/BioJournals.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Lesko02/BioJournals.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Lesko02/BioJournals.jl)


Simulators tracking mutations in a cell population, e.g., to study the
evolution of a solid tumor, must manipulate many DNA/RNA sequences at
a time.  It is obviously unfeasible to duplicate a full sequence every
time a mutation happens.  Any such simulator must therefore deal rely
on sophisticated data structures just tracking the *delta*'s from one
sequence to the next, while keeping them in a *tree-like* (or more
complex) data structure.


## Goal

The goal of the the `BioJournals` project is to provide a **Julia**
library that implements a form of *journaled sequence*, alongside a
proper API for its manipulation.  An important *key* identifying such
sequences will be *time*.

Similar data structures that have been proposed to solve similar
problems (cf., [SeqAn](https://www.seqan.de/data-structures/), but
they target different problems that the one at hand.  Plus, we would
like a **Julia** implementation to complement our other efforts in the
(solid tumors) simulation arena.

