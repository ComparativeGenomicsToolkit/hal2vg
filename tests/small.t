#!/usr/bin/env bash

BASH_TAP_ROOT=${VGDIR}/deps/bash-tap
. ${VGDIR}/deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH
PATH=../deps/hal:$PATH

plan tests 5

maf2hal small/small.maf small.hal
hal2vg small.hal > small.vg
vg view -j small.vg | jq . > small.json

is $(vg validate small.vg | wc -l) 0 "output vg validates"

is $(jq --argfile a small.json --argfile b small/truth.json -n '($a == $b)') true "output graph identical to manually verified truth graph"

hal2vg small.hal --protoChunk 500 > small_chunk500.vg
vg view -j small_chunk500.vg | jq . > small_chunk500.json

is $(jq --argfile a small_chunk500.json --argfile b small/truth.json -n '($a == $b)') true "output graph using --protoChunk 500 identical to manually verified truth graph"

rm -f small.vg small.json small_chunk500.vg small_chunk500.json

hal2vg small.hal --refGenome cat > small_cat.vg
vg view -j small_cat.vg | jq . > small_cat.json

is $(vg validate small_cat.vg | wc -l) 0 "output cat-referenced vg validates"

is $(jq --argfile a small_cat.json --argfile b small/truth_cat.json -n '($a == $b)') true "output cat-referenced graph identical to manually verified truth graph"

rm -f small_cat.vg small_cat.json

rm -f small.hal 
