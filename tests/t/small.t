#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=..:$PATH
PATH=../deps/hal/bin:$PATH

plan tests 9

maf2hal small/small.maf small.hal
hal2vg small.hal > small.vg
vg view -j small.vg | jq . > small.json

is $(vg validate small.vg | wc -l) 0 "output vg validates"

# jq craziness from https://stackoverflow.com/questions/31930041/using-jq-or-alternative-command-line-tools-to-compare-json-files
is $(jq --argjson a "$(cat small.json)" --argjson b "$(cat small/truth.json)" -n 'def post_recurse(f): def r: (f | select(. != null) | r), .; r; def post_recurse: post_recurse(.[]?); ($a | (post_recurse | arrays) |= sort) as $a | ($b | (post_recurse | arrays) |= sort) as $b | $a == $b') true "output graph identical to manually verified truth graph"

# --legacyConversion splits the pinch graph instead of planning the splits.  here the two happen to
# make the same graph, so it has to match the truth graph too
hal2vg small.hal --legacyConversion > small.legacy.vg
vg view -j small.legacy.vg | jq . > small.legacy.json
is $(jq --argjson a "$(cat small.legacy.json)" --argjson b "$(cat small/truth.json)" -n 'def post_recurse(f): def r: (f | select(. != null) | r), .; r; def post_recurse: post_recurse(.[]?); ($a | (post_recurse | arrays) |= sort) as $a | ($b | (post_recurse | arrays) |= sort) as $b | $a == $b') true "--legacyConversion output graph identical to manually verified truth graph"

rm -f small.vg small.json small.legacy.vg small.legacy.json

# a block with a member on the reverse strand, and up to four different bases at a position (the
# lower case t is the same base as T), in which the two conversions store some nodes in opposite
# orientations and can leave different nodes in two.  compare what neither changes: the spelled
# paths, and once unchopped, the edge count and each node's sequence on whichever strand sorts first
printf '##maf version=1\n\na score=0\n%s\n%s\n%s\n%s\n\na score=0\n%s\n%s\n%s\n' \
       's human.1  0 8 + 12 ACGTACGT' 's chimp.2  0 8 + 12 ACGAACGN' 's cat.3    2 8 - 10 ACCTACGt' \
       's dog.4    0 8 + 12 ACGCACGG' 's human.1  8 4 + 12 GATT' 's chimp.2  8 4 + 12 GACT' \
       's dog.4    8 2 + 12 GA--' > rev.maf
maf2hal rev.maf rev.hal > /dev/null
hal2vg rev.hal > rev.vg
hal2vg rev.hal --legacyConversion > rev.legacy.vg
unchopped() {
    vg mod -u "$1" | vg view - | awk -F'\t' '
        $1 == "S" { s = toupper($3); r = ""
                    for (i = length(s); i > 0; --i) { c = substr(s, i, 1); r = r (c == "A" ? "T" : c == "C" ? "G" : c == "G" ? "C" : c == "T" ? "A" : c) }
                    print "S", (s < r ? s : r) }
        $1 == "L" { ++edges }
        END { print "edges", edges }' | sort
}
is $(vg validate rev.vg | wc -l) 0 "a reversed member's output vg validates"
is "$(vg paths -F -x rev.vg)" "$(vg paths -F -x rev.legacy.vg)" "a reversed member's paths are spelled as --legacyConversion spells them"
is "$(unchopped rev.vg)" "$(unchopped rev.legacy.vg)" "a reversed member's unchopped graph is --legacyConversion's, up to orientation"

rm -f rev.maf rev.hal rev.vg rev.legacy.vg

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
