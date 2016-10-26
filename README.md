# hal2vg
Prototype code for converting [HAL](https://github.com/glennhickey/hal) to [vg](https://github.com/vgteam/vg).

(c) 2016 Glenn Hickey. See [LICENSE](https://github.com/glennhickey/hal2sg/blob/development/LICENSE) for details.

See also:
* [hal2sg](https://github.com/glennhickey/hal2sg): Convert  [HAL](https://github.com/glennhickey/hal) (output by [Cactus](https://github.com/glennhickey/progressiveCactus) and [CAMEL](https://github.com/adamnovak/sequence-graphs)) to [Side Graph SQL](https://github.com/ga4gh/schemas/wiki/Human-Genome-Variation-Reference-(HGVR)-Pilot-Project#graph-format)
* [sg2vg](https://github.com/glennhickey/sg2vg): Convert [Global Alliance (Side Graph) Server](https://github.com/ga4gh/schemas/wiki/Human-Genome-Variation-Reference-(HGVR)-Pilot-Project#graph-format) to [VG](https://github.com/ekg/vg)
* [vg2sg](https://github.com/glennhickey/vg2sg): Convert  [VG](https://github.com/ekg/vg) to [Side Graph SQL](https://github.com/ga4gh/schemas/wiki/Human-Genome-Variation-Reference-(HGVR)-Pilot-Project#graph-format)

## To do:
* More testing (though sanity checks built into hal2sg and sg2vg passing on lots of data so far).
* Repackage as part of Cacuts or vg?
* Output protobuf directly?
* Profile on very large graphs

## Algorithm

This tool is a composition of `hal2sg` and `sg2vg`.  It converts HAL into an in-memory version of the GA4GH Side Graph format, then exports that to vg.  The main difference is that it couts out the need to go through the SQL graph server format.  All dependencies are included as submodules. 

## Instructions

**Cloning:** Don't forget to clone submodules:

     git clone https://github.com/glennhickey/hal2vg.git --recursive

To run the converter:

	  hal2vg input.hal > output.json

or to go directly to vg:

     hal2vg input.hal | vg view -Jv > output.vg

To see all the options, run with no args or use `--help`.
