binPath=${rootPath}
libPath=${rootPath}

sonLibRootPath=deps/sonLib
sonLibPath=${sonLibRootPath}/lib

halRootPath=deps/hal
halPath=${halRootPath}/lib
halIncPath=${halRootPath}/api/inc

libbdsgPath=${rootPath}/deps/libbdsg-easy

jemallocPath=${rootPath}/deps/jemalloc

include  ${sonLibRootPath}/include.mk

# The pinch graph is millions of small, same-sized allocations.  glibc malloc puts an
# 8 byte header on each and rounds up to 16, so a 48 byte pinch segment really costs 64;
# jemalloc keeps its metadata out of band and has a 48 byte size class.  That is worth
# well over 15% of peak memory on a large conversion.
# Disabled on Mac (as in cactus) and under the address sanitizer, which brings its own
# allocator.  Set DISABLE_JEMALLOC=1 to link the system allocator instead.
# note: sonLib's include.mk sets SYS from uname
jemalloc = on
ifeq (${SYS},Darwin)
	jemalloc = off
endif
ifeq (${CGL_DEBUG},ultra)
	jemalloc = off
endif
ifdef DISABLE_JEMALLOC
	jemalloc = off
endif

CFLAGS += -I ${sonLibPath}  -I ${halPath} -I ${halIncPath}
CXXFLAGS += -std=c++14 -I ${sonLibPath}  -I ${halPath} -I ${halIncPath} -I ${libbdsgPath}/include -UNDEBUG
basicLibs = ${halPath}/libHal.a ${sonLibPath}/stPinchesAndCacti.a ${sonLibPath}/sonLib.a ${sonLibPath}/cuTest.a ${libbdsgPath}/lib/libbdsg.a ${libbdsgPath}/lib/libhandlegraph.a ${libbdsgPath}/lib/libsdsl.a ${libbdsgPath}/lib/libdivsufsort.a ${libbdsgPath}/lib/libdivsufsort64.a
ifeq (${jemalloc},on)
	basicLibs += ${jemallocPath}/lib/libjemalloc.a
	CXXFLAGS += -DHAVE_JEMALLOC
endif
# note the := : basicLibsDependencies is used as make prerequisites, so it must
# capture only the archives, before the -l flags below are appended
basicLibsDependencies := ${basicLibs}

# hal's hdf5Filters.cpp implements custom HDF5 lz4/zstd codecs, so libHal.a now
# carries LZ4_/ZSTD_ references.  h5c++'s wrapper config doesn't pull these in,
# so every binary linking libHal.a needs them spelled out.
basicLibs += -llz4 -lzstd

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

