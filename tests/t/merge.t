#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=../bin:$PATH
PATH=../deps/hal:$PATH

plan tests 10

maf2hal small/small.maf small.hal
maf2hal small/small2.maf small2.hal
halMergeChroms small.hal,small2.hal merged1.hal
halValidate merged1.hal
is $? 0 "halMergeChroms produces valid hal"
hal2fasta small.hal chimp > chimp.fa
hal2fasta small2.hal chimp >> chimp.fa
hal2fasta merged1.hal chimp > chimp.merge.fa
diff chimp.fa chimp.merge.fa
is $? 0 "halMergeChroms preserves chimp sequence"
hal2fasta small.hal cat > cat.fa
hal2fasta merged1.hal cat > cat.merge.fa
diff cat.fa cat.merge.fa
is $? 0 "halMergeChroms preserves cat sequence"
hal2vg small.hal | vg mod -O - | vg ids -s - > small.vg
hal2vg small2.hal | vg mod -O - | vg ids -s - > small2.vg
hal2vg merged1.hal | vg mod -O - | vg ids -s - > merged1.vg
vg view small.vg | sort > small.gfa
vg view small2.vg | sort > small2.gfa
vg find -x merged1.vg -p cat#3:1 -c 1000 | vg ids -s - | vg view - | sort | sed -e 's/_0//g' | sed -e 's/_1//g' | sed -e "s/human chimp cat/human cat chimp/g" > merged1.comp1.gfa
vg find -x merged1.vg -p cow#3:1 -c 1000 | vg ids -s - | vg view - | sort | sed -e 's/_0//g' | sed -e 's/_1//g' | sed -e "s/human cow chimp/chimp human cow/g" > merged1.comp2.gfa
diff small.gfa merged1.comp1.gfa
is $? 0 "First component of merged graph identical to first input graph"
diff small2.gfa merged1.comp2.gfa
is $? 0 "Second component of merged graph identical to second input graph"

rm -f small.hal small2.halsmall.vg small2.vg small.gfa small2.gfa
rm -f merged1.hal merged1.vg merged1.comp1.gfa merged1.comp2.gfa
rm -f chimp.fa chimp.merge.fa
rm -f cat.fa cat.merge.fa

### copy paste above but change order ###

maf2hal small/small.maf small.hal
maf2hal small/small2.maf small2.hal
halMergeChroms small2.hal,small.hal merged1.hal
halValidate merged1.hal
is $? 0 "halMergeChroms produces valid hal"
hal2fasta small2.hal chimp > chimp.fa
hal2fasta small.hal chimp >> chimp.fa
hal2fasta merged1.hal chimp > chimp.merge.fa
diff chimp.fa chimp.merge.fa
is $? 0 "halMergeChroms preserves chimp sequence"
hal2fasta small.hal cat > cat.fa
hal2fasta merged1.hal cat > cat.merge.fa
diff cat.fa cat.merge.fa
is $? 0 "halMergeChroms preserves cat sequence"
hal2vg small.hal | vg mod -O - | vg ids -s - > small.vg
hal2vg small2.hal | vg mod -O - | vg ids -s - > small2.vg
hal2vg merged1.hal | vg mod -O - | vg ids -s - > merged1.vg
vg view small.vg | sort > small.gfa
vg view small2.vg | sort > small2.gfa
vg find -x merged1.vg -p cat#3:1 -c 1000 | vg ids -s - | vg view - | sort | sed -e 's/_0//g' | sed -e 's/_1//g' | sed -e "s/human chimp cat/human cat chimp/g" > merged1.comp1.gfa
vg find -x merged1.vg -p cow#3:1 -c 1000 | vg ids -s - | vg view - | sort | sed -e 's/_0//g' | sed -e 's/_1//g' | sed -e "s/human cow chimp/chimp human cow/g" > merged1.comp2.gfa
diff small.gfa merged1.comp1.gfa
is $? 0 "First component of merged graph identical to first input graph"
diff small2.gfa merged1.comp2.gfa
is $? 0 "Second component of merged graph identical to second input graph"

rm -f small.hal small2.halsmall.vg small2.vg small.gfa small2.gfa
rm -f merged1.hal merged1.vg merged1.comp1.gfa merged1.comp2.gfa
rm -f chimp.fa chimp.merge.fa
rm -f cat.fa cat.merge.fa
