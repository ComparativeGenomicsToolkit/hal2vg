# hal2vg
[![Build Status](https://travis-ci.org/ComparativeGenomicsToolkit/hal2vg.svg?branch=master)](https://travis-ci.org/ComparativeGenomicsToolkit/hal2vg)

Convert [HAL](https://github.com/glennhickey/hal) to [vg](https://github.com/vgteam/vg)-compatible sequence graph.

Supports the three sequence graph formats in [libbdsg](https://github.com/vgteam/libbdsg):
* PackedGraph (default)
* ODGI
* HashGraph

## Algorithm

1. Each sequence in the HAL is added as a thread to a [Pinch Graph](https://github.com/ComparativeGenomicsToolkit/pinchesAndCacti).
2. Exact pairwise alignment blocks (no gaps or substitutions) are extracted from each branch in the HAL tree and "pinched" in the graph
3. For each branch, bases in the child that have substitutions in the parent (snps) are aligned across the tree using the column iterator and all exact matches are extracted and pinched.
4. Pinch graph is cleaned up by merging trivial joins
5. Each HAL sequence is traced through the pinch graph, adding nodes and edges to the output sequence graph.  A table is maintained to map each pinch graph block to a sequence graph node.

## Suggested Postprocessing:

*  Sort the output with `vg ids --sort`.  

## Installation

### Binary Release

You can download a standalone binary for the latest release [here](https://github.com/ComparativeGenomicsToolkit/hal2vg/releases).

### Compile From Source

You can use the the [Dockerfile](Dockerfile) as a guide to see how all dependencies are installed with `apt` on Ubuntu.  More details on installing HDF5 can be found in the [HAL README](https://github.com/ComparativeGenomicsToolkit/hal)

**Cloning:** Don't forget to clone submodules with the `--recursive` option:

     git clone https://github.com/glennhickey/hal2vg.git --recursive

**Compiling:**

     make

hal2vg links against [jemalloc](https://github.com/jemalloc/jemalloc), which is vendored as
a submodule and needs `autoconf` to build.  It saves a substantial amount of memory on large
conversions, since the pinch graph hal2vg builds internally is dominated by small, uniformly
sized allocations that jemalloc packs much more tightly than the system allocator.  It is
skipped automatically on macOS; build with `DISABLE_JEMALLOC=1 make` to use the system
allocator anywhere else.

## Usage

It is required to use the `--inMemory` option for all but trivial inputs.

`vg` has been tuned to work best on graphs with nodes chopped to at most 32 bases.  It is therefore recommended to use the `--chop 32` option.

```
hal2vg input.hal --inMemory --chop 32 --progress > output.pg
```

**Note**: The output graph is only readable by vg version 1.24.0 and greater.

### How blocks that disagree are converted

hal2vg pinches each genome to its parent along whole aligned segments, mismatches included, so a
block of the pinch graph can hold more than one base at a position.  By default, hal2vg works out
where each such block divides and how, and makes the nodes from that plan as it builds the graph,
without changing the pinch graph.  DNA is read in windows to do this: `--baseBuffer` sets their
total size in bytes (default 67108864), and is not a limit on memory use.

`--legacyConversion` instead splits those blocks in the pinch graph itself, as hal2vg used to,
which costs a pinch segment for every member of every piece and so far more memory on large
alignments.  The two make the same graph, apart from node IDs and a few nodes that the default
leaves in two where `--legacyConversion` makes one (`vg mod -u` removes the difference).  So
compare path metadata, spelled sequences and base-level homology rather than files or node
counts.  `make test-low-memory` does that for generated HALs with substitutions, indels,
reversals, paralogy, masked bases and path naming variants, and needs `vg` and Python 3.  For
complete PackedGraph outputs, see [tests/graph_equivalence.md](tests/graph_equivalence.md).

Copyright (C) 2020 by UCSC Computational Genomics Lab
