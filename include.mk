binPath=${rootPath}
libPath=${rootPath}

sonLibRootPath=deps/sonLib
sonLibPath=${sonLibRootPath}/lib

halRootPath=deps/hal
halPath=${halRootPath}/lib

hal2sgPath=${rootPath}/deps/hal2sg
sg2vgPath=${rootPath}/deps/sg2vg
rapidJsonPath=${sg2vgPath}/rapidjson
sgExportPath=${hal2sgPath}/sgExport

include  ${sonLibRootPath}/include.mk

cflags += -I ${sonLibPath}  -I ${halPath} -I ${sgExportPath} -I ${hal2sgPath} -I ${sg2vgPath} -I ${rapidJsonPath}/include
cppflags += -I ${sonLibPath}  -I ${halPath} -I ${sgExportPath} -I ${hal2sgPath} -I ${sg2vgPath} -I ${rapidJsonPath}/include
basicLibs = ${hal2sgPath}/libhal2sg.a ${sg2vgPath}/libsg2vg.a ${sgExportPath}/sgExport.a ${halPath}/halLiftover.a ${halPath}/halLib.a ${sonLibPath}/sonLib.a ${sonLibPath}/cuTest.a 
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


