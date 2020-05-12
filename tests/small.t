#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=../bin:$PATH
PATH=../deps/hal:$PATH

plan tests 4

maf2hal small/small.maf small.hal
hal2vg small.hal > small.vg
vg view -j small.vg | jq . > small.json

is $(vg validate small.vg | wc -l) 0 "output vg validates"

# jq craziness from https://stackoverflow.com/questions/31930041/using-jq-or-alternative-command-line-tools-to-compare-json-files
is $(jq --argfile a small.json --argfile b small/truth.json -n 'def post_recurse(f): def r: (f | select(. != null) | r), .; r; def post_recurse: post_recurse(.[]?); ($a | (post_recurse | arrays) |= sort) as $a | ($b | (post_recurse | arrays) |= sort) as $b | $a == $b') true "output graph identical to manually verified truth graph"

rm -f small.vg small.json

hal2vg small.hal --refGenome cat > small_cat.vg
vg view -j small_cat.vg | jq . > small_cat.json

is $(vg validate small_cat.vg | wc -l) 0 "output cat-referenced vg validates"

is $(jq --argfile a small_cat.json --argfile b small/truth_cat.json -n 'def post_recurse(f): def r: (f | select(. != null) | r), .; r; def post_recurse: post_recurse(.[]?); ($a | (post_recurse | arrays) |= sort) as $a | ($b | (post_recurse | arrays) |= sort) as $b | $a == $b') true "output cat-referenced graph identical to manually verified truth graph"

rm -f small_cat.vg small_cat.json

rm -f small.hal 
