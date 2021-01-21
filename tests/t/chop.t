#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=../bin:$PATH
PATH=../deps/hal:$PATH

plan tests 16

vg convert -g chop/tiny-flat.gfa -p > tiny-flat.vg
printf "x\t0\t100\n" > all.bed
chop-vg-paths tiny-flat.vg all.bed | vg view - | grep -v ^H > chopped-all.gfa
is "$(cat chopped-all.gfa | wc -l)" 0 "chopping everything clears out the graph"

rm -f all.bed chopped-all.gfa

printf "x\t0\t1\n" > ends.bed
printf "x\t48\t50\n" >> ends.bed
chop-vg-paths tiny-flat.vg ends.bed > chopped-ends.vg
is "$(vg paths -Ev chopped-ends.vg)" "x[1:48]	47" "chopping ends gives subpath in the middle with correct length"
is "$(vg stats -l chopped-ends.vg | awk '{print $2}')" "47" "chopping ends leaves correct number of bases"

rm -f ends.bed chopped-ends.vg

printf "x\t20\t25\n" > bits.bed
printf "x\t1\t5\n" >> bits.bed
printf "x\t10\t20\n" >> bits.bed
printf "x\t40\t49\n" >> bits.bed
chop-vg-paths tiny-flat.vg bits.bed > chopped-bits.vg
vg paths -Ev chopped-bits.vg | sed -e 's/\t/./g' >  bits.paths
is "$(cat bits.paths | wc -l)" "4" "correct number of paths obtained after merging consectuive intervals"
is "$(grep 'x\[0:1\].1' bits.paths | wc -l)" "1" "first bit found"
is "$(grep 'x\[5:10\].5' bits.paths | wc -l)" "1" "next bit found"
is "$(grep 'x\[25:40\].15' bits.paths | wc -l)" "1" "next bit after found"
is "$(grep 'x\[49:50\].1' bits.paths | wc -l)" "1" "last bit found"

rm -f bits.bed chopped-bits.vg bits.paths

rm -f tiny-flat.vg

########## flip path and repeat ##########

#vg convert -g chop/tiny-rev.gfa -p > tiny-rev.vg
vg convert -g chop/tiny-rev.gfa -o > tiny-rev.vg
printf "x\t0\t100\n" > all.bed
chop-vg-paths tiny-rev.vg all.bed | vg view - | grep -v ^H > chopped-all.gfa
is "$(cat chopped-all.gfa | wc -l)" 0 "chopping everything clears out the graph"

rm -f all.bed chopped-all.gfa

printf "x\t0\t1\n" > ends.bed
printf "x\t48\t50\n" >> ends.bed
chop-vg-paths tiny-rev.vg ends.bed > chopped-ends.vg
is "$(vg paths -Ev chopped-ends.vg)" "x[1:48]	47" "chopping ends gives subpath in the middle with correct length"
is "$(vg stats -l chopped-ends.vg | awk '{print $2}')" "47" "chopping ends leaves correct number of bases"

rm -f ends.bed chopped-ends.vg

printf "x\t20\t25\n" > bits.bed
printf "x\t1\t5\n" >> bits.bed
printf "x\t10\t20\n" >> bits.bed
printf "x\t40\t49\n" >> bits.bed
chop-vg-paths tiny-rev.vg bits.bed > chopped-bits.vg
vg paths -Ev chopped-bits.vg | sed -e 's/\t/./g' >  bits.paths
is "$(cat bits.paths | wc -l)" "4" "correct number of paths obtained after merging consectuive intervals"
is "$(grep 'x\[0:1\].1' bits.paths | wc -l)" "1" "first bit found"
is "$(grep 'x\[5:10\].5' bits.paths | wc -l)" "1" "next bit found"
is "$(grep 'x\[25:40\].15' bits.paths | wc -l)" "1" "next bit after found"
is "$(grep 'x\[49:50\].1' bits.paths | wc -l)" "1" "last bit found"

rm -f bits.bed chopped-bits.vg bits.paths

rm -f tiny-rev.vg
