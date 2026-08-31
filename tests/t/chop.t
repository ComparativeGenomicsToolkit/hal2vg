#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=..:$PATH
PATH=../deps/hal/bin:$PATH

plan tests 81

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
    # name the graph in the failure string: two invalid graphs must not compare equal, or an
    # assertion that both sides are unchanged passes when both sides in fact crashed
    if ! vg validate "$1" > /dev/null 2>&1; then echo "invalid-graph:$1"; return; fi
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

# --- -o/--out-bed ------------------------------------------------------------------------------

# -o reports what was clipped by diffing the input paths against the output ones, so it has to name
# the path a subpath came from and give coordinates in that path's space.  Two things were wrong,
# and neither shows on a graph whose paths carry no subranges -- where the key has no bracket to
# strip and the offset is 0 either way -- so this clips twice to get a subranged graph first.
{
    printf 'H\tVN:Z:1.0\n'
    for i in $(seq 12); do printf 'S\t%s\tACGTACGTAC\n' "$i"; done
    for i in $(seq 11); do printf 'L\t%s\t+\t%s\t+\t0M\n' "$i" "$((i+1))"; done
    seg="$(seq 12 | sed 's/$/+/' | paste -sd,)"
    printf 'P\tREF#0#chr1\t%s\t*\n' "$seg"
    printf 'P\tSAMP#0#chr1\t%s\t*\n' "$seg"
    printf 'P\t_MINIGRAPH_#0#a\t1+,2+\t*\n'
    printf 'P\t_MINIGRAPH_#0#b\t11+,12+\t*\n'
} > outbed.gfa
vg convert -g outbed.gfa -p > outbed.vg

# nodes 3-10 are off the minigraph, so -u 50 clips the middle and leaves SAMP as two subpaths
# carrying subranges: [0-20] and [100-120]
clip-vg outbed.vg -f -e REF -a _MINIGRAPH_ -u 50 > outbed-sub.vg
is "$(vg paths -L -x outbed-sub.vg 2>/dev/null | grep -c '^SAMP.*\[')" "2" "the first clip leaves subpaths carrying subranges"

# now clip those away and report it.  -m clips whole paths, so the two intervals are exactly the
# two subranges, and the BED must say so.
clip-vg outbed-sub.vg -f -e REF -m 25 -o outbed.bed > outbed-clipped.vg

# get_path_name() spells the subrange out as "name[start-end]", which files the record under a key
# no input path has and no downstream consumer can join on
is "$(grep -c '\[' outbed.bed)" "0" "-o names the path a subpath came from, not the subpath"

# the offset came from subrange.second -- where the subpath ends -- rather than .first, which put
# every interval a whole path-length to the right of the truth (20 40 / 120 140 here)
is "$(awk '$1 ~ /^SAMP/ {print $2"-"$3}' outbed.bed | sort -t- -k1,1n | tr '\n' ' ' | sed 's/ *$//')" \
   "0-20 100-120" "-o reports subpath intervals in the coordinates of the path they came from"

# the help has always advertised --no-orphan-filter; only --no-orphan_filter was ever wired up
clip-vg outbed.vg -f -e REF -a _MINIGRAPH_ -u 50 --no-orphan-filter > /dev/null 2>&1
is "$?" "0" "the --no-orphan-filter spelling the help documents is accepted"
clip-vg outbed.vg -f -e REF -a _MINIGRAPH_ -u 50 --no-orphan_filter > /dev/null 2>&1
is "$?" "0" "and the --no-orphan_filter spelling that shipped still works"

rm -f outbed.gfa outbed.vg outbed-sub.vg outbed.bed outbed-clipped.vg
# --- -k/--flank ------------------------------------------------------------------------------

# Build a synthetic graph from a spec of <count>:<nodelen>:<aligned 0|1> segments.  REF and SAMP
# run over every node; _MINIGRAPH_ covers only the aligned segments, so -a marks the rest unaligned.
mkgraph() { # $1=out.gfa  $2=spec
    awk -v spec="$2" 'BEGIN{
        print "H\tVN:Z:1.0"
        nseg=split(spec, S, " "); nid=0; allp=""
        for (s=1; s<=nseg; s++) {
            split(S[s], F, ":")
            for (i=0; i<F[1]; i++) {
                nid++; seq=""
                for (j=0; j<F[2]; j++) seq = seq substr("ACGT", (j%4)+1, 1)
                print "S\t" nid "\t" seq
                allp = allp (allp=="" ? "" : ",") nid "+"
                if (F[3]==1) mg[s] = (s in mg) ? mg[s] "," nid "+" : nid "+"
            }
        }
        for (i=1; i<nid; i++) print "L\t" i "\t+\t" (i+1) "\t+\t0M"
        print "P\tREF#0#chr1\t" allp "\t*"
        print "P\tSAMP#0#chr1\t" allp "\t*"
        for (s=1; s<=nseg; s++) if (s in mg) print "P\t_MINIGRAPH_#0#mg" s "\t" mg[s] "\t*"
    }' > "$1"
}

# the surviving SAMP fragments, as a sorted list of start-end subranges
samp_frags() { # $1=graph
    # name the graph, so two crashed runs never compare equal (see max_node_id)
    if ! vg validate "$1" > /dev/null 2>&1; then echo "invalid-graph:$1"; return; fi
    vg paths -E -x "$1" 2>/dev/null | awk '$1 ~ /^SAMP/ {print $1}' \
        | sed 's/.*\[//; s/\]//' | sort -t- -k1,1n | tr '\n' ' ' | sed 's/ *$//'
}

# aligned(200) unaligned(500) aligned-island(200) unaligned(500) aligned(200)
mkgraph island.gfa "2:100:1 5:100:0 2:100:1 5:100:0 2:100:1"
vg convert -g island.gfa -p > island.vg

clip-vg island.vg -f -e REF -a _MINIGRAPH_ -u 400 > island-u.vg
is "$(samp_frags island-u.vg)" "0-200 700-900 1400-1600" "an aligned island between two clipped runs survives -u"

# -k stops at that island rather than crossing it: the aligned nodes score negative and no
# unaligned sequence resumes beyond them before the next clipped interval.
clip-vg island.vg -f -e REF -a _MINIGRAPH_ -u 400 -k 1000 -T 0.25 > island-k.vg
is "$(samp_frags island-k.vg)" "0-200 700-900 1400-1600" "-k does not cross an aligned island"

rm -f island.gfa island.vg island-u.vg island-k.vg

# Two clipped runs separated by a gap that is mostly unaligned: aligned(20) unaligned(300)
# aligned(20).  Mott's rule walks each interval across its near aligned stretch because unaligned
# sequence resumes beyond it, so both extend into the gap and into each other.  chop_path cannot
# take overlapping intervals, so they have to be collapsed.
mkgraph collide.gfa "1:100:1 5:100:0 2:10:1 30:10:0 2:10:1 5:100:0 1:100:1"
vg convert -g collide.gfa -p > collide.vg

clip-vg collide.vg -f -e REF -a _MINIGRAPH_ -u 400 > collide-u.vg
is "$(samp_frags collide-u.vg)" "0-100 600-940 1440-1540" "without -k the sub-threshold gap survives"

clip-vg collide.vg -f -e REF -a _MINIGRAPH_ -u 400 -k 1000 -T 0.25 > collide-k.vg
is "$(samp_frags collide-k.vg)" "0-100 1440-1540" "-k extends two intervals into each other and they collapse"

rm -f collide.gfa collide.vg collide-u.vg collide-k.vg

# clean(400) fringe(390) anchor(10) unaligned(500) anchor(10) fringe(390) clean(400).  The fringe
# is unaligned too, but in a run of 390 -- under -u 400, so -u leaves it and only -k can reach it.
mkgraph flank.gfa "2:200:1 39:10:0 1:10:1 5:100:0 1:10:1 39:10:0 2:200:1"
vg convert -g flank.gfa -p > flank.vg

clip-vg flank.vg -f -e REF -a _MINIGRAPH_ -u 400 > flank-u.vg
is "$(samp_frags flank-u.vg)" "0-800 1300-2100" "without -k the sub-threshold fringe survives"

clip-vg flank.vg -f -e REF -a _MINIGRAPH_ -u 400 -k 1000 -T 0.25 > flank-k.vg
is "$(samp_frags flank-k.vg)" "0-400 1700-2100" "-k trims the unaligned fringe and stops at aligned sequence"
is "$(vg paths -E -x flank-k.vg 2>/dev/null | awk '$1 ~ /^REF/ {print $2}')" "2100" "-k leaves the reference intact"
is "$(max_node_id flank-u.vg)" "$(max_node_id flank-k.vg)" "-k mints no new node ids"

clip-vg flank.vg -f -e REF -a _MINIGRAPH_ -u 400 -k 100 -T 0.25 > flank-cap.vg
is "$(samp_frags flank-cap.vg)" "0-700 1400-2100" "-k never trims further than its cap"

rm -f flank.gfa flank.vg flank-u.vg flank-k.vg flank-cap.vg

# A flank only 20% unaligned by base (20bp unaligned against 80bp anchored).  A wholly unaligned
# flank goes at any threshold below 1, so it takes a mixed one to show that the threshold, rather
# than the mere presence of unaligned sequence, is what decides how far the trim runs.
mixed=""
for i in $(seq 10); do mixed="$mixed 1:20:0 1:80:1"; done
mkgraph mixed.gfa "2:200:1$mixed 5:100:0$mixed 2:200:1"
vg convert -g mixed.gfa -p > mixed.vg

# Below the threshold the whole flank goes.  The right edge stops at 2820 rather than 2900 because
# the trim ends on the last *unaligned* node -- running out over the trailing 80bp one lowers the
# score -- so the periodic pattern leaves the two sides a phase apart.
clip-vg mixed.vg -f -e REF -a _MINIGRAPH_ -u 400 -k 2000 -T 0.15 > mixed-lo.vg
is "$(samp_frags mixed-lo.vg)" "0-400 2820-3300" "-k trims a 20% unaligned flank when the threshold is below it"

# Above it neither side moves at all -- the first node each walk sees is 80bp of anchored sequence,
# which at this threshold outscores the 20bp of unaligned beyond it.  Assert that against -u alone
# rather than against a literal, so this stays a statement about -k doing nothing.  (The clipped
# run is 520bp, not 500: the 5:100:0 block is followed immediately by the next 1:20:0.)
clip-vg mixed.vg -f -e REF -a _MINIGRAPH_ -u 400 > mixed-none.vg
clip-vg mixed.vg -f -e REF -a _MINIGRAPH_ -u 400 -k 2000 -T 0.25 > mixed-hi.vg
is "$(samp_frags mixed-hi.vg)" "$(samp_frags mixed-none.vg)" "-k trims nothing when the threshold is above the flank's density"

# Everything above passes an explicit -T, so the default -- calibrating against the graph -- has
# never run.  Here unaligned sequence is spread evenly and the graph is only 3300bp long, so no base sits
# far enough from the clipped interval to measure a background rate from.  With nothing to compare
# against, calibration has to refuse: treating an unsampled background as a clean one puts the
# threshold at half the measured rate, which is below that rate by construction, so the walk never
# turns negative and the trim runs to its -k cap over ordinary sequence.
clip-vg mixed.vg -f -e REF -a _MINIGRAPH_ -u 400 -k 2000 > mixed-auto.vg 2> mixed-auto.err
is "$(samp_frags mixed-auto.vg)" "0-1400 1920-3300" "the default -T refuses to trim when it cannot sample a background"
is "$(grep -c 'too little background' mixed-auto.err)" "1" "and says so without needing -p"

# -0.0 is not less than zero, so an unnormalised "negative means calibrate" test would let it
# through to threshold 0 -- where a long node scores 0 rather than negative and the walk can never
# die.  Python's str(-0.0) is "-0.0", so a config can hand us one.
clip-vg mixed.vg -f -e REF -a _MINIGRAPH_ -u 400 -k 2000 -T -0 > mixed-negzero.vg 2>/dev/null
is "$(samp_frags mixed-negzero.vg)" "0-1400 1920-3300" "-T -0 calibrates like any other negative, rather than trimming flat out"

rm -f mixed.gfa mixed.vg mixed-lo.vg mixed-hi.vg mixed-none.vg mixed-auto.vg mixed-auto.err mixed-negzero.vg

# The same calibration where there is a background to measure: 300kb of clean 200bp nodes on each
# side, so sequence beyond the far window really is sampled, and the only unaligned sequence in
# the graph is the 390bp fringe flanking the clipped interval.  Runs on defaults -- no -T.
mkgraph calib.gfa "1500:200:1 39:10:0 1:10:1 5:100:0 1:10:1 39:10:0 1500:200:1"
vg convert -g calib.gfa -p > calib.vg

clip-vg calib.vg -f -e REF -a _MINIGRAPH_ -u 400 -k 5000 > calib-k.vg 2> calib-k.err
is "$(samp_frags calib-k.vg)" "0-300000 301300-601300" "the default -T trims exactly the fringe when it can sample a background"
is "$(grep -c 'using threshold 0.039' calib-k.err)" "1" "and reports the threshold it calibrated to"

rm -f calib.gfa calib.vg calib-k.vg calib-k.err

# -b intervals, unlike the -u ones every case above uses, need not begin or end on a node boundary.
# This graph is a palindrome and the interval is centred on it, both endpoints sitting 50bp inside
# a clean 1000bp node, so whatever -k trims off one side it must trim off the other.  Scoring the
# part of the straddled node that lies outside the interval is what keeps the two sides in step;
# charging the whole node on one side and skipping it entirely on the other does not.
mkgraph straddle.gfa "1:1000:1 5:20:0 1:1000:1 1:500:0 1:1000:1 5:20:0 1:1000:1"
vg convert -g straddle.gfa -p > straddle.vg
straddle_path="$(vg paths -E -x straddle.vg 2>/dev/null | awk '$1 ~ /^SAMP/ {print $1}')"
printf '%s\t2050\t2650\n' "$straddle_path" > straddle.bed

clip-vg straddle.vg -f -e REF -a _MINIGRAPH_ -b straddle.bed -k 3000 -T 0.2 > straddle-hi.vg
is "$(samp_frags straddle-hi.vg)" "0-2050 2650-4700" "-k charges the straddled node on both sides, so neither side trims"

clip-vg straddle.vg -f -e REF -a _MINIGRAPH_ -b straddle.bed -k 3000 -T 0.05 > straddle-lo.vg
is "$(samp_frags straddle-lo.vg)" "0-1000 3700-4700" "-k still trims both flanks, symmetrically, at a lower threshold"

rm -f straddle.gfa straddle.vg straddle.bed straddle-hi.vg straddle-lo.vg

# Two clipped intervals with a clean, fully anchored island between them.  The island is 0%
# unaligned, so -k has to stop at it: the sequence past it looks unaligned only because it is
# itself being clipped, and a walk allowed to score that would drag the boundary across it.
mkgraph nbr.gfa "2:200:1 5:100:0 3:100:1 25:20:0 2:200:1"
vg convert -g nbr.gfa -p > nbr.vg

clip-vg nbr.vg -f -e REF -a _MINIGRAPH_ -u 400 -k 600 -T 0.15 > nbr-k.vg
is "$(samp_frags nbr-k.vg)" "0-400 900-1200 1700-2100" "-k stops at a neighbouring clipped interval instead of scoring through it"

rm -f nbr.gfa nbr.vg nbr-k.vg

# Book-ended intervals are disjoint under half-open coordinates and load_bed accepts them.  The
# merge pass is there to collapse the overlaps -k can create, so it has to leave these alone:
# merging them would drop a breakpoint and quietly renumber what a plain -b run has always emitted.
mkgraph abut.gfa "10:100:1"
vg convert -g abut.gfa -p > abut.vg
abut_path="$(vg paths -E -x abut.vg 2>/dev/null | awk '$1 ~ /^SAMP/ {print $1}')"
printf '%s\t150\t250\n%s\t250\t350\n' "$abut_path" "$abut_path" > abut.bed

# -k is what puts the merge pass in play; without it the pass is skipped and this asserts nothing.
# Every node here is anchored, so no walk extends and the pass sees the two intervals as loaded.
clip-vg abut.vg -f -e REF -b abut.bed -k 1000 -T 0.25 > abut-b.vg
is "$(max_node_id abut-b.vg)" "13" "book-ended BED intervals keep the breakpoint they share"

rm -f abut.gfa abut.vg abut.bed abut-b.vg

# load_bed only rejects strict overlaps, so a file can carry an empty, inverted or off-the-end
# record.  Both walks would repair one into a plausible interval that chop_path then honours,
# turning a typo into deleted sequence -- and for a wholly off-path record, an abort into silence.
# -k has to leave anything it cannot reason about exactly as it arrived.
mkgraph degen.gfa "4:100:1 39:10:0 4:100:1"
vg convert -g degen.gfa -p > degen.vg
degen_path="$(vg paths -E -x degen.vg 2>/dev/null | awk '$1 ~ /^SAMP/ {print $1}')"

for rec in "600 600" "700 650" "1190 1300"; do
    printf '%s\t%s\t%s\n' "$degen_path" ${rec} > degen.bed
    clip-vg degen.vg -f -e REF -b degen.bed > degen-b.vg 2>/dev/null
    clip-vg degen.vg -f -e REF -a _MINIGRAPH_ -b degen.bed -k 1000 -T 0.25 > degen-k.vg 2>/dev/null
    is "$(samp_frags degen-k.vg)" "$(samp_frags degen-b.vg)" "-k leaves a [$rec] BED record alone"
done

# wholly off the path: master aborts, and -k must not turn that into a clean exit
printf '%s\t1500\t1600\n' "$degen_path" > degen.bed
clip-vg degen.vg -f -e REF -b degen.bed > /dev/null 2>&1; degen_rc=$?
clip-vg degen.vg -f -e REF -a _MINIGRAPH_ -b degen.bed -k 1000 -T 0.25 > /dev/null 2>&1
is "$?" "$degen_rc" "-k does not convert an off-path BED record into a clean exit"

rm -f degen.gfa degen.vg degen.bed degen-b.vg degen-k.vg

# --- option handling ---------------------------------------------------------------------------

mkgraph opt.gfa "10:100:1"
vg convert -g opt.gfa -p > opt.vg

# every node scores negative at or above 1, so no flank could ever be trimmed: say so up front
# rather than run to completion having done nothing
clip-vg opt.vg -f -e REF -a _MINIGRAPH_ -u 400 -k 1000 -T 1 > /dev/null 2> opt-t1.err
is "$?" "1" "-T at or above 1 is rejected instead of silently never trimming"

# NaN passes both signbit and >= 1, then makes every score NaN so no comparison against the
# running maximum is ever true: -k does nothing, exit 0, empty log.  Reject it instead.
clip-vg opt.vg -f -e REF -a _MINIGRAPH_ -u 400 -k 1000 -T nan > /dev/null 2> opt-nan.err
is "$?" "1" "-T nan is rejected rather than silently disabling -k"

# 0 is in the documented range and is a coherent request -- extend unconditionally -- but at it an
# aligned base scores 0 rather than negative, so the walk cannot die and every trim runs to the cap
clip-vg opt.vg -f -e REF -a _MINIGRAPH_ -u 400 -k 1000 -T 0 > /dev/null 2> opt-t0.err
is "$(grep -c 'gates nothing' opt-t0.err)" "1" "-T 0 warns that it gates nothing"

# -k needs intervals to act on, but a wrapper may pass it unconditionally and disable with 0,
# so this warns rather than failing a job that is otherwise fine
clip-vg opt.vg -f -e REF -k 1000 > /dev/null 2> opt-nosrc.err
is "$?" "0" "-k without an interval source is not fatal"
is "$(grep -c 'needs clipped intervals' opt-nosrc.err)" "1" "but it does warn"

rm -f opt.gfa opt.vg opt-t1.err opt-nan.err opt-t0.err opt-nosrc.err
