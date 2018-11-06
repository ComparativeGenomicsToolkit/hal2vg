# hal2vg
Prototype code for converting [HAL](https://github.com/glennhickey/hal) to [vg](https://github.com/vgteam/vg).

(c) 2016 Glenn Hickey. See [LICENSE](https://github.com/glennhickey/hal2vg/blob/master/LICENSE) for details.

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

#### vg

* [vg](https://github.com/vgteam/vg) must be downloaded and built before hal2vg
* Edit hal2vg/include.mk and make sure that VGDIR points to the correct vg directory

#### HDF5 1.10.1 with C++ API enabled

* Use build from [Progressive Cactus](https://github.com/glennhickey/progressiveCactus) by downloading and building Progressive Cactus *then* running `. environment` in the progressive cactus directory before building hal2vg. 

* Or Local install from source into DIR (do not need root password)  

     `mkdir DIR/hdf5`  
     `wget http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.gz`  
     `tar xzf  hdf5-1.10.1.tar.gz`  
     `cd hdf5-1.10.1`  
     `./configure --enable-cxx --prefix DIR/hdf5`  
     `make && make install` 
    
     Before building HAL, update the following environment variables:  
   
     `export PATH=DIR/hdf5/bin:${PATH}`  
     `export h5prefix=-prefix=DIR/hdf5`  

* Or From [Source](http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/):

     `wget http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.gz`  
     `tar xzf  hdf5-1.10.1.tar.gz`  
     `cd hdf5-1.10.1`  
     `./configure --enable-cxx`  
     `make && make install`  

* Or Using [MacPorts](http://www.macports.org/):   

    sudo port install hdf5 @1.10.1 +cxx

## Instructions

**Cloning:** Don't forget to clone submodules with the `--recursive` option:

     git clone https://github.com/glennhickey/hal2vg.git --recursive

**Setting your VG path:** Edit `include.mk` so that `VGDIR` points to where you've built [vg](https://github.com/vgteam/vg).  By default it will be `../vg`

**Compiling:**

     make

To run the converter:

	  hal2vg input.hal > output.vg

To see all the options, run with no args or use `--help`.

Note: The output vg may have nodes with sequence length up to 1MB, and will need to be chopped (ex `vg mod -X 32`) before indexing with `vg index`.   
