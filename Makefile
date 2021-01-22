# simplest possible to start.  dangerous since all header deps done manually. ne
rootPath = ./
include ./include.mk

all : hal2vg chop-vg-paths

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

cleanFast : 
	rm -f hal2vg hal2vg.o chop-vg-paths chop-vg-paths.o

clean :
	rm -f hal2vg hal2vg.o chop-vg-paths chop-vg-paths.o
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

chop-vg-paths.o : chop-vg-paths.cpp ${basicLibsDependencies}
	${cpp} ${CXXFLAGS} -I . chop-vg-paths.cpp -c

chop-vg-paths : chop-vg-paths.o ${basicLibsDependencies}
	${cpp} ${CXXFLAGS} -fopenmp -pthread chop-vg-paths.o  ${basicLibs}  -o chop-vg-paths

test :
	make
	cd tests && prove -v t
