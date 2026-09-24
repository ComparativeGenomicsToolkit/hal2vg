#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=..:$PATH

plan tests 10

# A 40kb reference of 40 1kb nodes.  ctgA is one contig whose first 3kb maps to the reference's
# end before the rest maps to its start: a 37kb deletion asserted by a single contig, with a 3kb
# cheaper side.  ctgB/C/D each map nodes 1-10 and then 25-34: the same 16kb deletion asserted by
# three contigs, cheaper side 10kb.  -m 0.01 budgets 370bp and 160bp respectively, so on its own
# it leaves every one of them in.  Sequences never get aligned here, so they are all the same.
seq=$(printf 'ACGT%.0s' $(seq 250))
{
    printf 'H\tVN:Z:1.0\n'
    for i in $(seq 1 40); do printf 'S\ts%d\t%s\tSN:Z:GRCh38#0#chrT\tSO:i:%d\tSR:i:0\n' $i "$seq" $(( (i-1)*1000 )); done
    for i in $(seq 1 39); do printf 'L\ts%d\t+\ts%d\t+\t0M\n' $i $((i+1)); done
} > filter-ref.gfa
paf() { # $1=contig $2=contig length $3=query start $4=node number
    printf 'id=%s\t%d\t%d\t%d\t+\tid=_MINIGRAPH_|s%d\t1000\t0\t1000\t1000\t1000\t60\ttp:A:P\n' "$1" "$2" "$3" $(( $3 + 1000 )) "$4"
}
{
    q=0; for n in 36 37 38; do paf HGA.1\|ctgA 32000 $q $n; q=$((q+1000)); done
    for n in $(seq 2 30); do paf HGA.1\|ctgA 32000 $q $n; q=$((q+1000)); done
    for c in HGB.1\|ctgB HGC.1\|ctgC HGD.1\|ctgD; do
        q=0; for n in $(seq 1 10); do paf $c 20000 $q $n; q=$((q+1000)); done
        for n in $(seq 25 34); do paf $c 20000 $q $n; q=$((q+1000)); done
    done
} > filter-test.paf
vg convert -r 0 -g filter-ref.gfa -p -T filter-ref.trans > filter-ref.vg

kept() { # $1=contig  $2..=options
    local c=$1; shift
    filter-paf-deletions filter-ref.vg filter-ref.trans filter-test.paf -d 5000 -m 0.01 "$@" 2> /dev/null | awk -v c="id=$c" '$1==c' | wc -l
}

is $(kept HGA.1\|ctgA) 32 "without -M the singleton deletion is left in"
is $(kept HGB.1\|ctgB) 20 "without -M the consensus deletion is left in"

is $(kept HGA.1\|ctgA -M 20000 -S 2) 29 "-M resolves the singleton by removing its 3kb side"
is $(kept HGB.1\|ctgB -M 20000 -S 2) 20 "-M leaves a deletion three contigs assert (support 3 >= 2)"
is $(kept HGC.1\|ctgC -M 20000 -S 2) 20 "...for every contig asserting it"
is $(kept HGD.1\|ctgD -M 20000 -S 2) 20 "...for every contig asserting it"

is $(kept HGB.1\|ctgB -M 20000 -S 4) 10 "-S above the support makes the consensus eligible too"
is $(kept HGA.1\|ctgA -M 20000 -S 4) 29 "...and the singleton still goes"

is $(kept HGA.1\|ctgA -M 2000 -S 2) 32 "-M below the cheaper side leaves the singleton alone"
is $(kept HGB.1\|ctgB -M 2000 -S 2) 20 "...and the consensus"

rm -f filter-ref.gfa filter-test.paf filter-ref.vg filter-ref.trans
