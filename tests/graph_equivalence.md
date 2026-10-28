# Comparing complete PackedGraph outputs

Compare outputs from the same HAL and the same selection/reference options,
with `--legacyConversion` as the only difference. Use raw hal2vg
outputs, before clipping, filtering or other graph transformations.

Build the optional verifier using the existing dependencies plus OpenSSL/libcrypto:

```bash
make -f Makefile -f verification.mk hal2vg-fingerprint
./hal2vg-fingerprint legacy.pg > legacy.jsonl
./hal2vg-fingerprint default.pg > default.jsonl
python3 tests/compare_graph_fingerprints.py legacy.jsonl default.jsonl
```

Run the two scans sequentially. Each scan loads one output graph plus node-level
anchors, a bounded-per-node sequence cache and edge coverage information. It does
not construct per-base graph objects or export full graph JSON/FASTA. Loading a
large PackedGraph still requires substantial memory.

The verifier compares complete case-sensitive path sequences and metadata, and
canonical oriented base-level homology maps. It also requires every path
transition to have an edge and every graph node/edge to be traversed by a path.
This last requirement scopes the verifier to path-covered raw exports; it is not
a general rejection criterion for graphs from other sources.

Anchoring shared bases to their first occurrence on lexically sorted paths and
coalescing consecutive mapping intervals makes the SHA256 fingerprints independent
of node numbering and equivalent linear subdivisions. Missing completion records
or any unequal path record cause the comparison to fail. Equivalence establishes
agreement between two conversions, not the absence of a defect shared by both.

Small positive and deliberately damaged graph fixtures test the verifier itself:

```bash
make -f Makefile -f verification.mk test-graph-fingerprints
```
