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
3. At most two HAL genomes are ever loaded at a time
4. Pinch graph is cleaned up by merging trivial joins
5. Each HAL sequence is traced through the pinch graph, adding nodes and edges to the output sequence graph.  A table is maintained to map each pinch graph block to a sequence graph node.

## Suggested Postprocessing:

*  Remove overlong nodes with `vg mod --chop 32`.  vg's mapping logic is mostly tuned and tested using a 32-base node length cap.
*  Sort the output with `vg ids --sort`.  

## Installation

You can use the the [Dockerfile](Dockerfile) as a guide to see how all dependencies are installed with `apt` on Ubuntu.  More details on installing HDF5 can be found in the [HAL README](https://github.com/ComparativeGenomicsToolkit/hal)

**Cloning:** Don't forget to clone submodules with the `--recursive` option:

     git clone https://github.com/glennhickey/hal2vg.git --recursive

**Compiling:**

     make

## Usage

It is required to use the `--inMemory` option for all but trivial inputs.

```
hal2vg input.hal --inMemory > output.pg
```

**Note**: The output graph is only readable by vg version 1.24.0 and greater.

`hal2vg` uses lots of memory to store both the entire pinch graph and sequence graph.  An alternative conversion pipeline can be found [here](https://github.com/ComparativeGenomicsToolkit/hal#pangenome-graph-export-gfa-and-vg).

Copyright (C) 2020 by UCSC Computational Genomics Lab

