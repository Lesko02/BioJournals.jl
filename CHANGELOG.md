# Changelog

## 1.0.0 — 2026-11-03

First stable API.

### Added

- Journals for `LongDNA{2}`, `LongDNA{4}`, `LongRNA{2}`, and `LongRNA{4}`.
- `sequence`, `sequences`, and `nsequences`. `sequence` extends `FASTX.sequence` and returns a `LongSequence`.
- `ReferenceInterval` and `locus`, using 1-based closed intervals. An empty interval anchors an insertion before reference base `first`.
- `read_journal` builds a `JournaledString` from one FASTA or FASTQ record. `write_fasta` writes derived sequences in index order.
- `edits_at` returns the edits of one derived sequence in time order.
- A single sequence index is accepted by `add_delta!`, as well as any integer vector.

### Changed

- Positions are 1-based coordinates on the reference. Replay walks edits in time order, so removing an edit does not require rebasing the others.
- Mutation payloads are typed (`Deletion`, `Insertion`, `Substitution`, `StructuralVariant`, `CopyNumber`).
- The journal is the source of truth. The tree exposed by `print_tree` is rebuilt from it.
- Allele sequences stored in a journal are copied. A later change to the caller's sequence does not rewrite history.
- A structural variant past the end of a 2-bit sequence is rejected when it is recorded. Those alphabets cannot encode a gap.
- FASTA and FASTQ loaders reject duplicate record identifiers.
- Public docstrings follow the Julia manual style, with checked examples.
- The README documents coordinates, alphabets, search, and file input.

### Removed

- `insert!` is no longer part of the package. It collided with `Base.insert!`. Sequence edits go through `add_delta!`.
- `build_tree!`, `delta_count`, `JSTNode`, `delete_at!`, `structure_variation!`, `copy_number_variation!`, and `approximate_findall` are no longer exported.
