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

hal2vg.o : hal2vg.cpp ${basicLibsDependencies}
	${cpp} ${CXXFLAGS} -I . hal2vg.cpp -c

${sonLibPath}/sonLib.a :
	cd deps/sonLib && make

${halPath}/libHal.a : ${sonLibPath}/sonLib.a
	cd deps/hal && make

${sonLibPath}/stPinchesAndCacti.a : ${sonLibPath}/sonLib.a
	cd deps/pinchesAndCacti && make

${libbdsgPath}/lib/libbdsg.a :
	cd deps/libbdsg-easy && make

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
