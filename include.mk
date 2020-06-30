binPath=${rootPath}
libPath=${rootPath}

sonLibRootPath=deps/sonLib
sonLibPath=${sonLibRootPath}/lib

halRootPath=deps/hal
halPath=${halRootPath}/lib
halIncPath=${halRootPath}/api/inc

libbdsgPath=${rootPath}/deps/libbdsg-easy

include  ${sonLibRootPath}/include.mk

CFLAGS += -I ${sonLibPath}  -I ${halPath} -I ${halIncPath}
CXXFLAGS += -std=c++11 -I ${sonLibPath}  -I ${halPath} -I ${halIncPath} -I ${libbdsgPath}/include -UNDEBUG
basicLibs = ${halPath}/libHal.a ${sonLibPath}/stPinchesAndCacti.a ${sonLibPath}/sonLib.a ${sonLibPath}/cuTest.a ${libbdsgPath}/lib/libbdsg.a ${libbdsgPath}/lib/libhandlegraph.a ${libbdsgPath}/lib/libsdsl.a ${libbdsgPath}/lib/libdivsufsort.a ${libbdsgPath}/lib/libdivsufsort64.a
basicLibsDependencies = ${basicLibs}

# hdf5 compilation is done through its wrappers.
# we can speficy our own (sonlib) compilers with these variables:
HDF5_CXX = ${cpp}
HDF5_CXXLINKER = ${cpp}
HDF5_CC = ${cxx}
HDF5_CCLINKER = ${cxx} 
cpp = h5c++ ${h5prefix}
cxx = h5cc ${h5prefix}

# add compiler flag and kent paths if udc is enabled
# relies on KENTSRC containing path to top level kent/ dir
# and MACHTYPE being specified
ifdef ENABLE_UDC
#  Find samtabix as in kent/src/inc/common.mk:
	ifeq (${SAMTABIXDIR},)
		SAMTABIXDIR = /hive/data/outside/samtabix/${MACHTYPE}
	endif

	basicLibs += ${KENTSRC}/src/lib/${MACHTYPE}/jkweb.a  ${SAMTABIXDIR}/libsamtabix.a -lssl -lcrypto
endif

