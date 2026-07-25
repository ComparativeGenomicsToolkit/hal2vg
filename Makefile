# simplest possible to start.  dangerous since all header deps done manually. ne
rootPath = ./
include ./include.mk

all : hal2vg clip-vg halRemoveDupes halMergeChroms halUnclip filter-paf-deletions count-vg-hap-cov

# Note: hdf5 from apt doesn't seem to work for static builds.  It should be installed
# from source and configured with "--enable-static --disable-shared", then have its
# bin put at the front of PATH
static:
	CFLAGS="$${CFLAGS} -static" \
	CXXFLAGS="$${CXXFLAGS} -static" \
	${MAKE} all

check-static: static
	if [ $(shell ls hal2vg clip-vg halRemoveDupes halMergeChroms halUnclip filter-paf-deletions count-vg-hap-cov | xargs ldd 2>& 1 | grep "not a dynamic" | wc -l) = $(shell ls hal2vg clip-vg halRemoveDupes halMergeChroms halUnclip filter-paf-deletions count-vg-hap-cov | wc -l) ] ; then\
		echo "ldd verified that all files in bin/ are static";\
	else\
		echo "ldd found dynamic linked binary in bin/";\
		exit 1;\
	fi

cleanFast : 
	rm -f hal2vg hal2vg.o clip-vg clip-vg.o halRemoveDupes halRemoveDupes.o halMergeChroms halMergeChroms.o halUnclip halUnclip.o filter-paf-deletions filter-paf-deletions.o count-vg-hap-cov.o count-vg-hap-cov

clean :
	rm -f hal2vg hal2vg.o clip-vg clip-vg.o halRemoveDupes halRemoveDupes.o halMergeChroms halMergeChroms.o halUnclip halUnclip.o filter-paf-deletions filter-paf-deletions.o
	cd deps/sonLib && make clean
	cd deps/pinchesAndCacti && make clean
	cd deps/hal && make clean
	cd deps/libbdsg-easy && make clean
	if [ -e deps/jemalloc/Makefile ] ; then cd deps/jemalloc && make clean ; fi

hal2vg.o : hal2vg.cpp ${basicLibsDependencies}
	${cpp} ${CXXFLAGS} -I . hal2vg.cpp -c

# These recurse into submodules whose own makefiles know how to rebuild themselves, but
# make still has to be told when to bother recursing.  Without the source lists below,
# editing a submodule and running make here silently relinks the stale archive.
sonLibSources = $(wildcard deps/sonLib/C/impl/*.c deps/sonLib/C/inc/*.h)
pinchSources = $(wildcard deps/pinchesAndCacti/impl/*.c deps/pinchesAndCacti/inc/*.h)
halSources = $(shell find deps/hal -name '*.cpp' -o -name '*.h' 2>/dev/null)

${sonLibPath}/sonLib.a : ${sonLibSources}
	cd deps/sonLib && make

${halPath}/libHal.a : ${sonLibPath}/sonLib.a ${halSources}
	cd deps/hal && make

${sonLibPath}/stPinchesAndCacti.a : ${sonLibPath}/sonLib.a ${pinchSources}
	cd deps/pinchesAndCacti && make

${libbdsgPath}/lib/libbdsg.a :
	cd deps/libbdsg-easy && make

# jemalloc is vendored (rather than taken from the system) so that the static release
# binaries are self contained and always get the same allocator: its size classes are
# what make the pinch graph cheap, so a different version could quietly change how much
# memory hal2vg needs.  --disable-libdl keeps it linkable into a fully static binary.
# jemalloc is configured with its own clean flags rather than the ones make static
# exports: its configure probes the number of significant virtual address bits by
# linking a test program, and that probe fails outright under -static.  The archive
# links into a static binary regardless, since -static is a link time flag.
${jemallocPath}/lib/libjemalloc.a :
	cd deps/jemalloc && \
	  env CFLAGS="-O3 -g" CXXFLAGS="-O3 -g" LDFLAGS="" LIBS="" ./autogen.sh && \
	  env CFLAGS="-O3 -g" CXXFLAGS="-O3 -g" LDFLAGS="" LIBS="" ./configure --disable-libdl && \
	  env CFLAGS="-O3 -g" CXXFLAGS="-O3 -g" LDFLAGS="" LIBS="" ${MAKE} build_lib_static

hal2vg : hal2vg.o ${basicLibsDependencies}
	${cpp} ${CXXFLAGS} -fopenmp -pthread hal2vg.o  ${basicLibs}  -o hal2vg

clip-vg.o : clip-vg.cpp ${basicLibsDependencies}
	${cpp} ${CXXFLAGS} -I . clip-vg.cpp -c

clip-vg : clip-vg.o ${basicLibsDependencies}
	${cpp} ${CXXFLAGS} -fopenmp -pthread clip-vg.o  ${basicLibs}  -o clip-vg

halRemoveDupes.o : halRemoveDupes.cpp ${basicLibsDependencies}
	${cpp} ${CXXFLAGS} -I . halRemoveDupes.cpp -c

halRemoveDupes : halRemoveDupes.o ${basicLibsDependencies}
	${cpp} ${CXXFLAGS} -fopenmp -pthread halRemoveDupes.o  ${basicLibs}  -o halRemoveDupes

halMergeChroms.o : halMergeChroms.cpp ${basicLibsDependencies}
	${cpp} ${CXXFLAGS} -I . halMergeChroms.cpp -c

halMergeChroms : halMergeChroms.o ${basicLibsDependencies}
	${cpp} ${CXXFLAGS} -fopenmp -pthread halMergeChroms.o  ${basicLibs}  -o halMergeChroms

halUnclip.o : halUnclip.cpp subpaths.h ${basicLibsDependencies}
	${cpp} ${CXXFLAGS} -I . halUnclip.cpp -c

halUnclip : halUnclip.o ${basicLibsDependencies} 
	${cpp} ${CXXFLAGS} -fopenmp -pthread halUnclip.o  ${basicLibs}  -o halUnclip

filter-paf-deletions.o : filter-paf-deletions.cpp subpaths.h paf.hpp ${basicLibsDependencies}
	${cpp} ${CXXFLAGS} -I . filter-paf-deletions.cpp -c

filter-paf-deletions : filter-paf-deletions.o ${basicLibsDependencies} 
	${cpp} ${CXXFLAGS} -fopenmp -pthread filter-paf-deletions.o  ${basicLibs}  -o filter-paf-deletions

count-vg-hap-cov.o : count-vg-hap-cov.cpp ${basicLibsDependencies}
	${cpp} ${CXXFLAGS} -I . count-vg-hap-cov.cpp -c

count-vg-hap-cov : count-vg-hap-cov.o ${basicLibsDependencies}
	${cpp} ${CXXFLAGS} -fopenmp -pthread count-vg-hap-cov.o  ${basicLibs}  -o count-vg-hap-cov

test :
	make
	cd tests && prove -v t
