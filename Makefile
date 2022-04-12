# simplest possible to start.  dangerous since all header deps done manually. ne
rootPath = ./
include ./include.mk

all : hal2vg clip-vg halRemoveDupes halMergeChroms halUnclip filter-paf-deletions

# Note: hdf5 from apt doesn't seem to work for static builds.  It should be installed
# from source and configured with "--enable-static --disable-shared", then have its
# bin put at the front of PATH
static:
	CFLAGS="$${CFLAGS} -static" \
	CXXFLAGS="$${CXXFLAGS} -static" \
	${MAKE} all

check-static: static
ifeq ($(shell ldd hal2vg | grep "not a dynamic" | wc -l), $(shell ls hal2vg | wc -l))
	$(info ldd verified that hal2vg static)
else
	$(error ldd found dnymaic linked dependency in hal2vg)
endif
ifeq ($(shell ldd clip-vg | grep "not a dynamic" | wc -l), $(shell ls clip-vg | wc -l))
	$(info ldd verified that clip-vg static)
else
	$(error ldd found dnymaic linked dependency in clip-vg)
endif
ifeq ($(shell ldd halRemoveDupes | grep "not a dynamic" | wc -l), $(shell ls halRemoveDupes | wc -l))
	$(info ldd verified that halRemoveDupes static)
else
	$(error ldd found dnymaic linked dependency in halRemoveDupes)
endif
ifeq ($(shell ldd halMergeChroms | grep "not a dynamic" | wc -l), $(shell ls halMergeChroms | wc -l))
	$(info ldd verified that halMergeChroms static)
else
	$(error ldd found dnymaic linked dependency in halMergeChroms)
endif
ifeq ($(shell ldd halUnclip | grep "not a dynamic" | wc -l), $(shell ls halUnclip | wc -l))
	$(info ldd verified that halUnclip static)
else
	$(error ldd found dnymaic linked dependency in halUnclip)
endif
ifeq ($(shell ldd filter-paf-deletions | grep "not a dynamic" | wc -l), $(shell ls filter-paf-deletions | wc -l))
	$(info ldd verified that filter-paf-deletions static)
else
	$(error ldd found dnymaic linked dependency in filter-paf-deletions)
endif


cleanFast : 
	rm -f hal2vg hal2vg.o clip-vg clip-vg.o halRemoveDupes halRemoveDupes.o halMergeChroms halMergeChroms.o halUnclip halUnclip.o filter-paf-deletions filter-paf-deletions.o

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

test :
	make
	cd tests && prove -v t
