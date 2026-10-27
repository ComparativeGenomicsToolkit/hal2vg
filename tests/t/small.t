#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=..:$PATH
PATH=../deps/hal/bin:$PATH

plan tests 7

maf2hal small/small.maf small.hal
hal2vg small.hal > small.vg
vg view -j small.vg | jq . > small.json

is $(vg validate small.vg | wc -l) 0 "output vg validates"

# jq craziness from https://stackoverflow.com/questions/31930041/using-jq-or-alternative-command-line-tools-to-compare-json-files
is $(jq --argjson a "$(cat small.json)" --argjson b "$(cat small/truth.json)" -n 'def post_recurse(f): def r: (f | select(. != null) | r), .; r; def post_recurse: post_recurse(.[]?); ($a | (post_recurse | arrays) |= sort) as $a | ($b | (post_recurse | arrays) |= sort) as $b | $a == $b') true "output graph identical to manually verified truth graph"

# --lowMemory never splits the pinch graph, so it can leave a node broken in two where the default
# conversion joins it.  once both are unchopped there is no such difference left
hal2vg small.hal --lowMemory > small.low.vg
is $(vg validate small.low.vg | wc -l) 0 "--lowMemory output vg validates"
is "$(vg mod -u small.low.vg | vg ids -s - | vg view - | sort)" "$(vg mod -u small.vg | vg ids -s - | vg view - | sort)" "--lowMemory output is the default output once both are unchopped"

rm -f small.vg small.json small.low.vg

rm -f small.hal 

# genome names are PanSN SAMPLE#HAPLOTYPE (or just SAMPLE, for haplotype 0), and a '.' is only a
# haplotype in the old SAMPLE.HAPLOTYPE spelling of a HAL made before cactus used them
maf2hal small/small.maf names.hal
printf 'human\tHG002#1\nchimp\tKOLF2.1J\ncat\tGCA_054883195.1#0\n' > names.tsv
halRenameGenomes names.hal names.tsv > /dev/null
hal2vg names.hal > names.vg
is "$(vg paths -L -v names.vg | sort | tr '\n' ' ')" "GCA_054883195.1#0#3 HG002#1#1 KOLF2.1J#0#2 " "PanSN genome names make PanSN paths"

printf 'HG002#1\tHG002.2\n' > names.tsv
halRenameGenomes names.hal names.tsv > /dev/null
hal2vg names.hal > names.vg
is "$(vg paths -L -v names.vg | sort | tr '\n' ' ')" "GCA_054883195.1#0#3 HG002#2#1 KOLF2.1J#0#2 " "the old SAMPLE.HAPLOTYPE spelling is still read"

printf 'HG002.2\tHG002#x\n' > names.tsv
halRenameGenomes names.hal names.tsv > /dev/null
hal2vg names.hal > names.vg 2> /dev/null
isnt $? 0 "a # that is not SAMPLE#HAPLOTYPE is an error"

rm -f names.hal names.tsv names.vg
