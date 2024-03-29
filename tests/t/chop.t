#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=../bin:$PATH
PATH=../deps/hal:$PATH

plan tests 18

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

rm -f tiny-fr.vg tiny-fr.fa tiny-fr-forwardized.vg tiny-fr-forwardized.fa tiny-fr-forwardized.fa
