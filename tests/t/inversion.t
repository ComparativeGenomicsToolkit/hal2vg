#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=..:$PATH

plan tests 33

# A 2kb reference of 20 100bp nodes.  HGX walks r1-r3 forward, a 50bp junction, r15 down to r5
# backwards (a 1.1kb inversion), 1.5kb of sequence nothing else aligns to, then r16-r20 forward:
# -u 1000 clips the 1.5kb and with it the inversion's right junction, leaving the left one.  HGY
# has the same inversion with a 50bp junction at each end.  HGZ is the reference.  clip-vg never
# aligns anything, so the sequences are all the same.
node() { printf 'S\t%s\t%s\n' "$1" "$(printf 'ACGT%.0s' $(seq $(( $2 / 4 + 1 ))) | head -c $2)"; }
{
    printf 'H\tVN:Z:1.0\n'
    for i in $(seq 1 20); do node r$i 100; done
    node a1 50; node a2 50; node a3 50; node u1 500; node u2 500; node u3 500
    for i in $(seq 1 19); do printf 'L\tr%d\t+\tr%d\t+\t0M\n' $i $((i+1)); done
    printf 'L\tr3\t+\ta1\t+\t0M\nL\ta1\t+\tr15\t-\t0M\nL\tr5\t-\tu1\t+\t0M\nL\tu1\t+\tu2\t+\t0M\nL\tu2\t+\tu3\t+\t0M\nL\tu3\t+\tr16\t+\t0M\n'
    printf 'L\tr3\t+\ta2\t+\t0M\nL\ta2\t+\tr15\t-\t0M\nL\tr5\t-\ta3\t+\t0M\nL\ta3\t+\tr16\t+\t0M\n'
    rev=$(for i in $(seq 15 -1 5); do printf 'r%d-,' $i; done)
    fwd=$(for i in $(seq 1 20); do printf 'r%d+,' $i; done)
    printf 'P\tGRCh38#0#chrT\t%s\t*\n' "${fwd%,}"
    printf 'P\tHGX#1#ctg\tr1+,r2+,r3+,a1+,%su1+,u2+,u3+,r16+,r17+,r18+,r19+,r20+\t*\n' "$rev"
    printf 'P\tHGY#1#ctg\tr1+,r2+,r3+,a2+,%sa3+,r16+,r17+,r18+,r19+,r20+\t*\n' "$rev"
    printf 'P\tHGZ#1#ctg\t%s\t*\n' "${fwd%,}"
} > inversion.gfa
vg convert -g inversion.gfa -p > inversion.vg

# subpaths of a sample after clipping, as "start-end:steps:reverse-steps" per line
fragments() { # $1=graph  $2=sample
    vg view "$1" | awk -F'\t' -v s="$2" '$1=="W" && $2==s {w=$7; n=gsub(/[<>]/,"",w); r=gsub(/</,"",$7); print $5"-"$6":"n":"r}' | sort -t- -k1,1n | tr '\n' ' '
}

clip-vg inversion.vg -e GRCh38 -u 1000 > noI.vg 2> /dev/null
is "$(fragments noI.vg HGX)" "0-1450:15:11 2950-3450:5:0 " "without -I the one-sided inversion keeps its surviving junction"
is "$(fragments noI.vg HGY)" "0-2000:21:11 " "without -I the two-sided inversion is one path"

clip-vg inversion.vg -e GRCh38 -u 1000 -I 500 > withI.vg 2> /dev/null
is "$(fragments withI.vg HGX)" "0-350:4:0 350-1450:11:11 2950-3450:5:0 " "-I severs the surviving junction: the inverted run is its own subpath and nothing is clipped"
is "$(fragments withI.vg HGY)" "0-2000:21:11 " "-I leaves a two-sided inversion alone"
is "$(fragments withI.vg HGZ)" "0-2000:20:0 " "-I leaves a forward path alone"
is "$(vg validate withI.vg 2>&1)" "graph: valid" "the graph is still valid"

clip-vg inversion.vg -e GRCh38 -u 1000 -I 2000 > big.vg 2> /dev/null
is "$(fragments big.vg HGX)" "0-1450:15:11 2950-3450:5:0 " "-I above the run's size leaves it"

clip-vg inversion.vg -u 1000 -I 500 > /dev/null 2> noref.err
isnt "$(grep -c 'requires -e' noref.err)" 0 "-I without -e is refused"

# A haplotype that carries a segment both forward and inverted.  HGR walks r1, 2kb nothing else
# aligns to, r3-r2 backwards (a 200bp inversion whose left junction that 2kb clips away), r4, r2
# forward again, r5, another such 2kb and a second such inversion, r7-r6, then r8-r9.  -I splits
# the path at each run's surviving junction and clips nothing, so no node is divided and the
# reference keeps its nine steps.  -f as cactus runs it.
{
    printf 'H\tVN:Z:1.0\n'
    for i in $(seq 1 9); do node r$i 100; done
    node x 2000; node y 2000
    for i in $(seq 1 8); do printf 'L\tr%d\t+\tr%d\t+\t0M\n' $i $((i+1)); done
    printf 'L\tr1\t+\tx\t+\t0M\nL\tx\t+\tr3\t-\t0M\nL\tr2\t-\tr4\t+\t0M\nL\tr4\t+\tr2\t+\t0M\nL\tr2\t+\tr5\t+\t0M\n'
    printf 'L\tr5\t+\ty\t+\t0M\nL\ty\t+\tr7\t-\t0M\nL\tr6\t-\tr8\t+\t0M\n'
    printf 'P\tGRCh38#0#chrT\tr1+,r2+,r3+,r4+,r5+,r6+,r7+,r8+,r9+\t*\n'
    printf 'P\tHGR#1#ctg\tr1+,x+,r3-,r2-,r4+,r2+,r5+,y+,r7-,r6-,r8+,r9+\t*\n'
} > revisit.gfa
vg convert -g revisit.gfa -p > revisit.vg

clip-vg revisit.vg -f -e GRCh38 -u 1000 -I 150 > revisit.out.vg 2> revisit.err
is "$?" 0 "-I on a path that carries a segment both ways"
is "$(grep -c 'Severed 2 of 2' revisit.err)" 1 "both one-sided inversions were severed"
is "$(fragments revisit.out.vg HGR)" "0-100:1:0 2100-2300:2:2 2300-2600:3:0 4600-4800:2:2 4800-5000:2:0 " "each run is split off at its junction with nothing clipped"
is "$(fragments revisit.out.vg GRCh38)" "0-900:9:0 " "no node was divided: the reference keeps its nine steps"
is "$(vg validate revisit.out.vg 2>&1)" "graph: valid" "the graph is still valid"

# The chopper on that same path, given mid-node cuts through -b: one base at the end of each
# inverted run (the junction bases the first version of -I clipped) plus the two 2kb stretches
# (vg names the P line path HGR#1#ctg#0).
# The first cut divides r2 during its reverse visit; the chopper used to take each step's length
# as it went, so the later forward visit of r2 counted only the piece that kept r2's id, every
# later offset was 99 bases short, the r6 cut missed its node and clip-vg died on an assertion.
printf 'HGR#1#ctg#0\t100\t2100\nHGR#1#ctg#0\t2299\t2300\nHGR#1#ctg#0\t2600\t4600\nHGR#1#ctg#0\t4799\t4800\n' > revisit.bed
clip-vg revisit.vg -f -e GRCh38 -b revisit.bed > revisit.bed.vg 2> /dev/null
is "$?" 0 "mid-node cuts on a path that revisits a divided node"
is "$(fragments revisit.bed.vg HGR)" "0-100:1:0 2100-2299:2:2 2300-2600:4:0 4600-4799:2:2 4800-5000:2:0 " "each cut landed in its own node and the revisit is intact"
is "$(fragments revisit.bed.vg GRCh38)" "0-900:11:0 " "the reference keeps the cut bases, as two extra steps"
is "$(vg validate revisit.bed.vg 2>&1)" "graph: valid" "the graph is still valid"

# Severing is decided against the graph, not the path alone.  HGS carries HGX's inversion with both
# junctions (through a1 and a3), so HGX's surviving junction edge a1->r15- is walked by a two-sided
# path: severing HGX would change no topology, only fragment it, and -I leaves it alone.  HGX2 is a
# second one-sided carrier (its own unaligned tail v1-v3 gets clipped): when every path walking the
# junction edge is one-sided there, all of them are severed and the edge goes.
{
    grep -v '^P' inversion.gfa
    printf 'P\tGRCh38#0#chrT\t%s\t*\n' "${fwd%,}"
    printf 'P\tHGX#1#ctg\tr1+,r2+,r3+,a1+,%su1+,u2+,u3+,r16+,r17+,r18+,r19+,r20+\t*\n' "$rev"
    printf 'P\tHGY#1#ctg\tr1+,r2+,r3+,a2+,%sa3+,r16+,r17+,r18+,r19+,r20+\t*\n' "$rev"
    printf 'P\tHGS#1#ctg\tr1+,r2+,r3+,a1+,%sa3+,r16+,r17+,r18+,r19+,r20+\t*\n' "$rev"
} > shared.gfa
vg convert -g shared.gfa -p > shared.vg
clip-vg shared.vg -e GRCh38 -u 1000 -I 500 > shared.out.vg 2> shared.err
is "$(grep -c 'Severed 0 of 3 .*(1 left alone: another path walks the junction edge)' shared.err)" 1 "a junction edge a two-sided path walks is left alone (three runs found: HGX, HGY and HGS)"
is "$(fragments shared.out.vg HGX)" "0-1450:15:11 2950-3450:5:0 " "...and the one-sided path is not fragmented for nothing"
is "$(fragments shared.out.vg HGS)" "0-2000:21:11 " "...while the two-sided path is untouched"
{
    grep -v '^P' inversion.gfa
    node v1 500; node v2 500; node v3 500
    printf 'L\tr5\t-\tv1\t+\t0M\nL\tv1\t+\tv2\t+\t0M\nL\tv2\t+\tv3\t+\t0M\nL\tv3\t+\tr16\t+\t0M\n'
    printf 'P\tGRCh38#0#chrT\t%s\t*\n' "${fwd%,}"
    printf 'P\tHGX#1#ctg\tr1+,r2+,r3+,a1+,%su1+,u2+,u3+,r16+,r17+,r18+,r19+,r20+\t*\n' "$rev"
    printf 'P\tHGX2#1#ctg\tr1+,r2+,r3+,a1+,%sv1+,v2+,v3+,r16+,r17+,r18+,r19+,r20+\t*\n' "$rev"
} > twice.gfa
vg convert -g twice.gfa -p > twice.vg
clip-vg twice.vg -e GRCh38 -u 1000 -I 500 > twice.out.vg 2> twice.err
is "$(grep -c 'Severed 2 of 2 reverse-strand runs >= 500 bp that clipping had left with one junction$' twice.err)" 1 "two one-sided carriers of the same junction edge are both severed"
is "$(fragments twice.out.vg HGX)" "0-350:4:0 350-1450:11:11 2950-3450:5:0 " "...the first"
is "$(fragments twice.out.vg HGX2)" "0-350:4:0 350-1450:11:11 2950-3450:5:0 " "...and the second"

# A 2kb reference of 20 100bp nodes.  HGV is a contig the assembler emitted on the reverse strand:
# it walks r20 down to r11 backwards, then r9,r10 forward (a 200bp inversion, both junctions intact),
# then r8 down to r1 backwards.  Its reverse flanks are not inversions: -I must not sever them.
# HGW is the same reverse contig whose last 600bp are forward and run to the contig end (a one-sided
# inversion from the contig's point of view): its forward run has one intact junction and is severed.
# HGN is a forward contig with a one-step reverse blip in its flank (r3) and a 600bp inversion at its
# end: the blip is noise, not a run, and must not make the flank look like a bounded island.
{
    printf 'H\tVN:Z:1.0\n'
    for i in $(seq 1 20); do node r$i 100; done
    for i in $(seq 1 19); do printf 'L\tr%d\t+\tr%d\t+\t0M\n' $i $((i+1)); done
    printf 'L\tr11\t-\tr9\t+\t0M\nL\tr10\t+\tr8\t-\t0M\nL\tr11\t-\tr5\t+\t0M\n'
    printf 'L\tr2\t+\tr3\t-\t0M\nL\tr3\t-\tr4\t+\t0M\nL\tr10\t+\tr20\t-\t0M\n'
    top=$(for i in $(seq 20 -1 11); do printf 'r%d-,' $i; done)
    bot=$(for i in $(seq 8 -1 1); do printf 'r%d-,' $i; done)
    printf 'P\tGRCh38#0#chrT\t%s\t*\n' "${fwd%,}"
    printf 'P\tHGV#1#ctg\t%sr9+,r10+,%s\t*\n' "$top" "${bot%,}"
    printf 'P\tHGW#1#ctg\t%sr5+,r6+,r7+,r8+,r9+,r10+\t*\n' "$top"
    printf 'P\tHGN#1#ctg\tr1+,r2+,r3-,r4+,r5+,r6+,r7+,r8+,r9+,r10+,r20-,r19-,r18-,r17-,r16-,r15-\t*\n'
} > revcontig.gfa
vg convert -g revcontig.gfa -p > revcontig.vg

clip-vg revcontig.vg -f -e GRCh38 -u 1000 -I 150 > revcontig.out.vg 2> revcontig.err
is "$?" 0 "-I on a reverse-strand contig"
is "$(grep -c 'Severed 2 of 4' revcontig.err)" 1 "of the four reverse runs, the two whose forward neighbour reaches a contig end are severed"
is "$(fragments revcontig.out.vg HGV)" "0-2000:20:18 " "a reverse contig with a 200bp forward island stays in one piece"
is "$(fragments revcontig.out.vg HGW)" "0-1000:10:10 1000-1600:6:0 " "a reverse contig whose forward run reaches its end is split at the run's intact junction"
is "$(fragments revcontig.out.vg HGN)" "0-1000:10:1 1000-1600:6:6 " "a one-step reverse blip in the flank does not stop a one-sided inversion from being severed"
is "$(fragments revcontig.out.vg GRCh38)" "0-2000:20:0 " "no node was divided"
is "$(vg validate revcontig.out.vg 2>&1)" "graph: valid" "the graph is still valid"

clip-vg revcontig.vg -f -e GRCh38 -u 1000 -I 250 > big.vg 2> big.err
is "$(grep -c 'Severed 2 of 4' big.err)" 1 "-I above the island's size changes nothing: the flanks are the runs and the island is what guards them"
is "$(fragments big.vg HGN)" "0-1000:10:1 1000-1600:6:6 " "and the one-sided inversion is still severed"
is "$(fragments big.vg HGV)" "0-2000:20:18 " "and leaves the contig alone"

rm -f inversion.gfa inversion.vg noI.vg withI.vg big.vg big.err noref.err revisit.gfa revisit.vg revisit.out.vg revisit.err revisit.bed revisit.bed.vg
rm -f shared.gfa shared.vg shared.out.vg shared.err twice.gfa twice.vg twice.out.vg twice.err revcontig.gfa revcontig.vg revcontig.out.vg revcontig.err
