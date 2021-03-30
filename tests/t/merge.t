#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=../bin:$PATH
PATH=../deps/hal:$PATH

plan tests 4

maf2hal small/small.maf small.hal
maf2hal small/small2.maf small2.hal
halMergeChroms small.hal,small2.hal merged1.hal
hal2vg small.hal --noAncestors | vg ids -s - > small.vg
hal2vg small2.hal --noAncestors | vg ids -s - > small2.vg
hal2vg merged1.hal --noAncestors | vg ids -s - > merged1.vg
vg view small.vg  > small.gfa
vg view small2.vg > small2.gfa
vg find -x merged1.vg -p cat.3:1 -c 1000 | vg view - > merged1.comp1.gfa
vg find -x merged1.vg -p cow.3:1 -c 1000 | vg view - > merged1.comp2.gfa
diff small.gfa merged1.comp1.gfa
is $? 0 "First component of merged graph identical to first input graph"
diff small2.gfa merged1.comp2.gfa
is $? 0 "Second component of merged graph identical to second input graph"

# repeat with merge in different order
halMergeChroms small2.hal,small.hal merged2.hal
hal2vg merged2.hal --noAncestors | vg ids -s - > merged2.vg
vg find -x merged2.vg -p cat.3:1 -c 1000 | vg view - > merged2.comp1.gfa
vg find -x merged2.vg -p cow.3:1 -c 1000 | vg view - > merged2.comp2.gfa
diff small.gfa merged1.comp1.gfa
is $? 0 "First component of merged graph2 identical to first input graph"
diff small2.gfa merged1.comp2.gfa
is $? 0 "Second component of merged graph2 identical to second input graph"

rm -f small.hal small2.halsmall.vg small2.vg small.gfa small2.gfa
rm -f merged1.hal merged1.vg merged1.comp1.gfa merged1.comp2.gfa
rm -f merged2.hal merged2.vg merged2.comp1.gfa merged2.comp2.gfa
