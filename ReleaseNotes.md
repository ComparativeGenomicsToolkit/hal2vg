# Release 1.0.1   2020-09-07

This release contains a bugfix required to use the subsetting options such as `--targetGenomes` and `--ignoreGenomes` without crashing.  

# Release 1.0.0   2020-08-19

This release uses Cactus's Pinch Graph library to create the sequence graph.

Notable Changes:
 - Bespoke STL-based structures and algorimths replaced with Pinch Graph library
 - SNPS are aligned through using the column iterator instead of tables
 - Much more performant than original implementaiton, but still only tested on smallish graphs
