# hal2vg
[![Build Status](https://travis-ci.org/ComparativeGenomicsToolkit/hal2vg.svg?branch=master)](https://travis-ci.org/ComparativeGenomicsToolkit/hal2vg)

Prototype code for converting [HAL](https://github.com/glennhickey/hal) to [vg](https://github.com/vgteam/vg).

See also:
* [hal2sg](https://github.com/glennhickey/hal2sg): Convert  [HAL](https://github.com/glennhickey/hal) (output by [Cactus](https://github.com/glennhickey/progressiveCactus) and [CAMEL](https://github.com/adamnovak/sequence-graphs)) to [Side Graph SQL](https://github.com/ga4gh/schemas/wiki/Human-Genome-Variation-Reference-(HGVR)-Pilot-Project#graph-format)
* [sg2vg](https://github.com/glennhickey/sg2vg): Convert [Global Alliance (Side Graph) Server](https://github.com/ga4gh/schemas/wiki/Human-Genome-Variation-Reference-(HGVR)-Pilot-Project#graph-format) to [VG](https://github.com/ekg/vg)
* [vg2sg](https://github.com/glennhickey/vg2sg): Convert  [VG](https://github.com/ekg/vg) to [Side Graph SQL](https://github.com/ga4gh/schemas/wiki/Human-Genome-Variation-Reference-(HGVR)-Pilot-Project#graph-format)

## To do:
* More testing (though sanity checks built into hal2sg and sg2vg passing on lots of data so far).
* Repackage as part of Progressive Cactus
* Profile on very large graphs

## Algorithm

This tool is a composition of `hal2sg` and `sg2vg`.  It converts HAL into an in-memory version of the GA4GH Side Graph format, then exports that to vg.  The main difference is that it removes the need to go through the SQL graph server format as an intermediate.  All dependencies except HDF5 and VG are included as submodules. 

## Installing Dependencies

#### HDF5 1.10.1 with C++ API enabled

* Using apt (Ubuntu 18.04)

    `sudo apt install libhdf5-dev`

* Using [MacPorts](http://www.macports.org/):   

    `sudo port install hdf5 @1.10.1 +cxx`

* From [Source](http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/):

     `wget http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.gz`  
	  `tar xzf  hdf5-1.10.1.tar.gz`  
     `cd hdf5-1.10.1`  
	  `./configure --enable-cxx`  
	  `make && make install`  

* Local install from source into DIR (do not need root password)  

     `mkdir DIR/hdf5`
     `wget http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.gz`
	  `tar xzf  hdf5-1.10.1.tar.gz`  
     `cd hdf5-1.10.1`  
     `./configure --enable-cxx --prefix DIR/hdf5`  
     `make && make install` 
    
     Before building HAL, update the following environment variables:  
   
     `export PATH=DIR/hdf5/bin:${PATH}`  
     `export h5prefix=-prefix=DIR/hdf5`  
 
     or set these in include.local.mk.

    If you are using older version of HDF5, such as installed on Centos,
    you may need to set 
    
    `export CXX_ABI_DEF='-D_GLIBCXX_USE_CXX11_ABI=1'
    
    If you get undefined functions base on string type with errors about
    `std::__cxx11::basic_string` vs `std::basic_string`.

## Instructions

**Cloning:** Don't forget to clone submodules with the `--recursive` option:

     git clone https://github.com/glennhickey/hal2vg.git --recursive

**Compiling:**

     make

To run the converter:

	  hal2vg input.hal > output.pg

To see all the options, run with no args or use `--help`.

Note: The output graph may have nodes with sequence length up to 1MB, and will need to be chopped (ex `vg mod -X 32`) before indexing with `vg index`.

Note: The output graph is only readable by vg version 1.24.0 and greater.

(c) 2016 Glenn Hickey. See [LICENSE](https://github.com/glennhickey/hal2vg/blob/master/LICENSE) for details.
