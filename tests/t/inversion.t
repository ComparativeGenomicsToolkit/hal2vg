#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=..:$PATH

plan tests 8

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
is "$(fragments withI.vg HGX)" "0-300:3:0 350-1450:11:11 2950-3450:5:0 " "-I severs the surviving junction: the inverted run is its own subpath"
is "$(fragments withI.vg HGY)" "0-2000:21:11 " "-I leaves a two-sided inversion alone"
is "$(fragments withI.vg HGZ)" "0-2000:20:0 " "-I leaves a forward path alone"
is "$(vg validate withI.vg 2>&1)" "graph: valid" "the graph is still valid"

clip-vg inversion.vg -e GRCh38 -u 1000 -I 2000 > big.vg 2> /dev/null
is "$(fragments big.vg HGX)" "0-1450:15:11 2950-3450:5:0 " "-I above the run's size leaves it"

clip-vg inversion.vg -u 1000 -I 500 > /dev/null 2> noref.err
isnt "$(grep -c 'requires -e' noref.err)" 0 "-I without -e is refused"

rm -f inversion.gfa inversion.vg noI.vg withI.vg big.vg noref.err
