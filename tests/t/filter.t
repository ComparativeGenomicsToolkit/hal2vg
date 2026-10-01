#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=..:$PATH

plan tests 12

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

# A deletion that -m resolves in a later pass must not be acted on by the support pass.  On a 50kb
# reference whose last node is 1.5kb, ctgS maps A=s21, then F=s50, then s1-s5: two big deletions,
# A->F (28kb, cheaper side A at 1kb) and F->s1 (48kb, cheaper side F at 1.5kb).  -m 0.034 budgets
# 952bp for the first, so it is left in and recorded, and 1632bp for the second, so F goes -- which
# also resolves the first, as A now sits next to s1-s5.  A stale record of A->F used to survive into
# the support pass and take A out too.
{
    printf 'H\tVN:Z:1.0\n'
    for i in $(seq 1 50); do l=1000; [ $i -eq 50 ] && l=1500; printf 'S\ts%d\t%s\tSN:Z:GRCh38#0#chrT\tSO:i:%d\tSR:i:0\n' $i "$(printf 'ACGT%.0s' $(seq $((l/4))))" $(( (i-1)*1000 )); done
    for i in $(seq 1 49); do printf 'L\ts%d\t+\ts%d\t+\t0M\n' $i $((i+1)); done
} > stale-ref.gfa
vg convert -r 0 -g stale-ref.gfa -p -T stale-ref.trans > stale-ref.vg
{
    paf HGS.1\|ctgS 7500 0 21
    printf 'id=HGS.1|ctgS\t7500\t1000\t2500\t+\tid=_MINIGRAPH_|s50\t1500\t0\t1500\t1500\t1500\t60\ttp:A:P\n'
    q=2500; for n in $(seq 1 5); do paf HGS.1\|ctgS 7500 $q $n; q=$((q+1000)); done
} > stale.paf
stale() { filter-paf-deletions stale-ref.vg stale-ref.trans stale.paf -d 20000 -m 0.034 "$@" 2> /dev/null | cut -f6 | tr '\n' ' '; }
is "$(stale)" "id=_MINIGRAPH_|s21 id=_MINIGRAPH_|s1 id=_MINIGRAPH_|s2 id=_MINIGRAPH_|s3 id=_MINIGRAPH_|s4 id=_MINIGRAPH_|s5 " "-m removes F and leaves A"
is "$(stale -M 5000 -S 2)" "$(stale)" "-M/-S changes nothing: the deletion A was recorded for was resolved by then"

rm -f filter-ref.gfa filter-test.paf filter-ref.vg filter-ref.trans stale-ref.gfa stale-ref.vg stale-ref.trans stale.paf
