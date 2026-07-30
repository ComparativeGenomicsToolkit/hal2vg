#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=..:$PATH
PATH=../deps/hal/bin:$PATH

plan tests 47

# how many steps of a given path still point backwards.  Reports instead of counting if the
# graph is unusable or the path is missing, so that a crashed run leaving an empty output
# can't be mistaken for a graph with no reverse steps.
ref_rev_steps() { # $1=graph  $2=path name
    if ! vg validate "$1" > /dev/null 2>&1; then echo "invalid-graph"; return; fi
    if [ "$(vg view "$1" | awk -v p="$2" '$1=="P" && $2==p' | wc -l)" -ne 1 ]; then
        echo "no-such-path"; return
    fi
    vg view "$1" | awk -v p="$2" '$1=="P" && $2==p {print $3}' | tr ',' '\n' | grep -c -- '-$'
}

# same, across every path in the graph
all_rev_steps() { # $1=graph
    if ! vg validate "$1" > /dev/null 2>&1; then echo "invalid-graph"; return; fi
    vg view "$1" | awk '$1=="P" {print $3}' | tr ',' '\n' | grep -c -- '-$'
}

# the largest node id in the graph.  clip-vg's second invocation in cactus builds the clipped
# graph and has to stay inside the id space the first one established, so it must never mint one.
max_node_id() { # $1=graph
    if ! vg validate "$1" > /dev/null 2>&1; then echo "invalid-graph"; return; fi
    vg view "$1" | awk '$1=="S" {if ($2+0 > m) m = $2+0} END {print m+0}'
}

# forwardizing must never alter what the paths spell out
same_path_seqs() { # $1=before  $2=after
    vg paths -Fv "$1" | sort > .seq-before.fa
    vg paths -Fv "$2" | sort > .seq-after.fa
    diff -q .seq-before.fa .seq-after.fa > /dev/null
    echo $?
    rm -f .seq-before.fa .seq-after.fa
}

vg convert -g chop/tiny-flat.gfa -p > tiny-flat.vg
printf "x\t0\t100\n" > all.bed
clip-vg tiny-flat.vg -b all.bed | vg view - | grep -v ^H > chopped-all.gfa
is "$(cat chopped-all.gfa | wc -l)" 0 "chopping everything clears out the graph"

rm -f all.bed chopped-all.gfa

printf "y\t0\t100\n" > none.bed
clip-vg tiny-flat.vg -b none.bed | vg view - | grep -v ^H > chopped-none.gfa
vg view tiny-flat.vg | grep -v ^H > orig.gfa
diff chopped-none.gfa orig.gfa
is "$?" 0 "chopping nothing doesn't change graph"

rm -f none.bed chopped-none.gfa orig.gfa

printf "x\t0\t1\n" > ends.bed
printf "x\t48\t50\n" >> ends.bed
clip-vg -n tiny-flat.vg -b ends.bed > chopped-ends.vg
is "$(vg paths -Ev chopped-ends.vg)" "x[1-48]	47" "chopping ends gives subpath in the middle with correct length"
is "$(vg stats -l chopped-ends.vg | awk '{print $2}')" "47" "chopping ends leaves correct number of bases"

rm -f ends.bed chopped-ends.vg

printf "x\t20\t25\n" > bits.bed
printf "x\t1\t5\n" >> bits.bed
printf "x\t10\t20\n" >> bits.bed
printf "x\t40\t49\n" >> bits.bed
clip-vg -n tiny-flat.vg -b bits.bed > chopped-bits.vg
vg paths -Ev chopped-bits.vg | sed -e 's/\t/./g' >  bits.paths
is "$(cat bits.paths | wc -l)" "4" "correct number of paths obtained after merging consectuive intervals"
is "$(grep 'x\[0-1\].1' bits.paths | wc -l)" "1" "first bit found"
is "$(grep 'x\[5-10\].5' bits.paths | wc -l)" "1" "next bit found"
is "$(grep 'x\[25-40\].15' bits.paths | wc -l)" "1" "next bit after found"
is "$(grep 'x\[49-50\].1' bits.paths | wc -l)" "1" "last bit found"

rm -f bits.bed chopped-bits.vg bits.paths

rm -f tiny-flat.vg

########## flip path and repeat ##########

vg convert -g chop/tiny-rev.gfa -p > tiny-rev.vg
#vg convert -g chop/tiny-rev.gfa -o > tiny-rev.vg
printf "x\t0\t100\n" > all.bed
clip-vg tiny-rev.vg -b all.bed | vg view - | grep -v ^H > chopped-all.gfa
is "$(cat chopped-all.gfa | wc -l)" 0 "chopping everything clears out the graph"

rm -f all.bed chopped-all.gfa

printf "x\t0\t1\n" > ends.bed
printf "x\t48\t50\n" >> ends.bed
clip-vg -n tiny-rev.vg -b ends.bed > chopped-ends.vg
is "$(vg paths -Ev chopped-ends.vg)" "x[1-48]	47" "chopping ends gives subpath in the middle with correct length"
is "$(vg stats -l chopped-ends.vg | awk '{print $2}')" "47" "chopping ends leaves correct number of bases"

rm -f ends.bed chopped-ends.vg

printf "x\t20\t25\n" > bits.bed
printf "x\t1\t5\n" >> bits.bed
printf "x\t10\t20\n" >> bits.bed
printf "x\t40\t49\n" >> bits.bed
clip-vg -n tiny-rev.vg -b bits.bed > chopped-bits.vg
vg paths -Ev chopped-bits.vg | sed -e 's/\t/./g' >  bits.paths
is "$(cat bits.paths | wc -l)" "4" "correct number of paths obtained after merging consectuive intervals"
is "$(grep 'x\[0-1\].1' bits.paths | wc -l)" "1" "first bit found"
is "$(grep 'x\[5-10\].5' bits.paths | wc -l)" "1" "next bit found"
is "$(grep 'x\[25-40\].15' bits.paths | wc -l)" "1" "next bit after found"
is "$(grep 'x\[49-50\].1' bits.paths | wc -l)" "1" "last bit found"

rm -f bits.bed chopped-bits.vg bits.paths

rm -f tiny-rev.vg

# quick test for forwardization
vg convert -g chop/tiny-fr.gfa -p > tiny-fr.vg
vg paths -Fv tiny-fr.vg > tiny-fr.fa
clip-vg tiny-fr.vg -e x -p > tiny-fr-forwardized.vg
vg paths -Fv tiny-fr-forwardized.vg > tiny-fr-forwardized.fa
diff tiny-fr.fa tiny-fr-forwardized.fa
is "$?" 0  "fowawrsization does not affect path sequence"

# the diff above only compares path sequence, which flipping preserves by construction.
# check that the reference actually came out forward.
is "$(ref_rev_steps tiny-fr-forwardized.vg x)" "0" "reference path x is forwardized"

rm -f tiny-fr.vg tiny-fr.fa tiny-fr-forwardized.vg tiny-fr-forwardized.fa

########## forwardizing a plain reverse reference node ##########

vg convert -g chop/ref-rev-simple.gfa -p > ref-rev.vg
clip-vg ref-rev.vg -e x > ref-rev-fwd.vg
is "$?" "0" "a reverse reference node forwardizes without error"
is "$(ref_rev_steps ref-rev-fwd.vg x)" "0" "reference path x is forwardized"
is "$(same_path_seqs ref-rev.vg ref-rev-fwd.vg)" "0" "forwardizing does not affect path sequence"

rm -f ref-rev.vg ref-rev-fwd.vg

########## reference cycles: -c must never abort ##########
# -c is used by minigraph-cactus when reference cycles are allowed.  It has to let
# non-forward reference nodes through rather than failing.

# every visit to the cycled node is backwards, so one flip forwardizes all of them
vg convert -g chop/ref-cycle-all-rev.gfa -p > cyc-rev.vg
clip-vg cyc-rev.vg -e x -c > cyc-rev-fwd.vg
is "$?" "0" "-c on an all-reverse reference cycle does not error"
is "$(ref_rev_steps cyc-rev-fwd.vg x)" "0" "an all-reverse reference cycle is forwardized"
is "$(same_path_seqs cyc-rev.vg cyc-rev-fwd.vg)" "0" "forwardizing a reference cycle preserves path sequence"

rm -f cyc-rev.vg cyc-rev-fwd.vg

# the cycled node is also the path's last step -- this used to walk the path traversal
# off its end and abort inside libbdsg
vg convert -g chop/ref-cycle-tail.gfa -p > cyc-tail.vg
clip-vg cyc-tail.vg -e x -c > cyc-tail-fwd.vg
is "$?" "0" "-c on a reference cycle ending on the cycled node does not crash"
is "$(ref_rev_steps cyc-tail-fwd.vg x)" "0" "a reference cycle on the last step is forwardized"
is "$(same_path_seqs cyc-tail.vg cyc-tail-fwd.vg)" "0" "forwardizing a tail cycle preserves path sequence"

rm -f cyc-tail.vg cyc-tail-fwd.vg

# the reference walks the node both ways round, so no orientation suits every visit.
# the reverse step has to survive into the output rather than being an error.
vg convert -g chop/ref-cycle-mixed.gfa -p > cyc-mix.vg
clip-vg cyc-mix.vg -e x -c > cyc-mix-fwd.vg
is "$?" "0" "-c on a mixed-orientation reference cycle does not error"
is "$(ref_rev_steps cyc-mix-fwd.vg x)" "1" "the un-forwardizable reference step passes through"
is "$(same_path_seqs cyc-mix.vg cyc-mix-fwd.vg)" "0" "a mixed reference cycle preserves path sequence"

# ... but without -c the same graph is still rejected
clip-vg cyc-mix.vg -e x > cyc-mix-nc.vg 2> cyc-mix-nc.err
is "$?" "1" "a reference cycle without -c is still an error"
is "$(grep -c 'Cycle detected' cyc-mix-nc.err)" "1" "the reference cycle error names the cause"

rm -f cyc-mix.vg cyc-mix-fwd.vg cyc-mix-nc.vg cyc-mix-nc.err

# two reference contigs sharing one node in opposite orientations -- not a cycle at all
vg convert -g chop/ref-two-contigs.gfa -p > two-ctg.vg
clip-vg two-ctg.vg -e x -c > two-ctg-fwd.vg
is "$?" "0" "-c on two reference contigs sharing a node does not error"
is "$(all_rev_steps two-ctg-fwd.vg)" "1" "the shared node stays reverse on one contig"

rm -f two-ctg.vg two-ctg-fwd.vg

########## non-reference cycles are always allowed ##########

vg convert -g chop/nonref-cycle.gfa -p > nr-cyc.vg
clip-vg nr-cyc.vg -e x > nr-cyc-fwd.vg
is "$?" "0" "an all-reverse non-reference cycle needs no -c"

vg convert -g chop/nonref-cycle-mixed.gfa -p > nr-cycm.vg
clip-vg nr-cycm.vg -e x > nr-cycm-fwd.vg
is "$?" "0" "a mixed-orientation non-reference cycle needs no -c"

rm -f nr-cyc.vg nr-cyc-fwd.vg nr-cycm.vg nr-cycm-fwd.vg

########## forwardizing nodes no path visits forward ##########
# a contig aligned in reverse leaves a run of nodes every path walks backwards.
# nothing downstream can flip them, so clip-vg does it here.

vg convert -g chop/nonref-rev.gfa -p > nr-rev.vg
clip-vg nr-rev.vg -e x -F > nr-rev-fwd.vg
is "$?" "0" "a reverse-aligned non-reference contig forwardizes without error"
is "$(ref_rev_steps nr-rev-fwd.vg samp.ctg)" "0" "nodes no path visits forward are flipped"
is "$(same_path_seqs nr-rev.vg nr-rev-fwd.vg)" "0" "flipping non-reference nodes preserves path sequence"

# -F mints node ids, so it must stay opt-in: cactus runs clip-vg a second time to build the
# clipped graph, and that run has to keep the id space of the first one
clip-vg nr-rev.vg -e x > nr-rev-noF.vg
is "$(ref_rev_steps nr-rev-noF.vg samp.ctg)" "2" "without -F, nodes no path visits forward are left alone"
is "$(vg stats -N nr-rev.vg)" "$(vg stats -N nr-rev-noF.vg)" "without -F the node count is unchanged"
is "$(max_node_id nr-rev.vg)" "$(max_node_id nr-rev-noF.vg)" "without -F no new node ids are minted"

clip-vg nr-rev.vg -F > nr-rev-badF.err 2>&1
is "$?" "1" "-F without -e is an error"

rm -f nr-rev-noF.vg nr-rev-badF.err

# ... and a prefix matching no path is an error, caught before anything is clipped
clip-vg nr-rev.vg -e no-such-prefix > nr-rev-none.vg 2> nr-rev-none.err
is "$?" "1" "a reference prefix matching no path is an error"
is "$(grep -c 'No path name begins with' nr-rev-none.err)" "1" "the error names the missing prefix"

# the case that mattered: with -u every base looks unaligned when there is no reference,
# and the orphan filter used to remove the whole graph and exit 0
clip-vg nr-rev.vg -e no-such-prefix -u 5 > nr-rev-u.vg 2> nr-rev-u.err
is "$?" "1" "a mistyped -e with -u fails instead of silently emitting an empty graph"

rm -f nr-rev.vg nr-rev-fwd.vg nr-rev-none.vg nr-rev-none.err nr-rev-u.vg nr-rev-u.err
