#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=..:$PATH
PATH=../deps/hal/bin:$PATH

plan tests 5

maf2hal small/small.maf small.hal
hal2vg small.hal > small.vg
vg view -j small.vg | jq . > small.json

is $(vg validate small.vg | wc -l) 0 "output vg validates"

# jq craziness from https://stackoverflow.com/questions/31930041/using-jq-or-alternative-command-line-tools-to-compare-json-files
is $(jq --argjson a "$(cat small.json)" --argjson b "$(cat small/truth.json)" -n 'def post_recurse(f): def r: (f | select(. != null) | r), .; r; def post_recurse: post_recurse(.[]?); ($a | (post_recurse | arrays) |= sort) as $a | ($b | (post_recurse | arrays) |= sort) as $b | $a == $b') true "output graph identical to manually verified truth graph"

rm -f small.vg small.json

rm -f small.hal 

# genome names are PanSN SAMPLE#HAPLOTYPE (or just SAMPLE, for haplotype 0), and a '.' is never a
# haplotype: cactus only reads the old SAMPLE.HAPLOTYPE spelling in its input seqfile
maf2hal small/small.maf names.hal
printf 'human\tHG002#1\nchimp\tKOLF2.1J\ncat\tGCA_054883195.1\n' > names.tsv
halRenameGenomes names.hal names.tsv > /dev/null
hal2vg names.hal > names.vg
is "$(vg paths -L -v names.vg | sort | tr '\n' ' ')" "GCA_054883195.1#0#3 HG002#1#1 KOLF2.1J#0#2 " "PanSN genome names make PanSN paths"

printf 'HG002#1\tHG002.2\nGCA_054883195.1\tGCA_054883195.1#0\n' > names.tsv
halRenameGenomes names.hal names.tsv > /dev/null
hal2vg names.hal > names.vg
is "$(vg paths -L -v names.vg | sort | tr '\n' ' ')" "GCA_054883195.1#0#3 HG002.2#0#1 KOLF2.1J#0#2 " "a '.' is part of the sample, and #0 is haplotype 0"

printf 'HG002.2\tHG002#x\n' > names.tsv
halRenameGenomes names.hal names.tsv > /dev/null
hal2vg names.hal > names.vg 2> /dev/null
isnt $? 0 "a # that is not SAMPLE#HAPLOTYPE is an error"

rm -f names.hal names.tsv names.vg
