# simplest possible to start.  dangerous since all header deps done manually. ne
rootPath = ./
include ./include.mk

all : hal2vg

cleanFast : 
	rm -f hal2vg hal2vg.o

clean :
	rm -f hal2vg hal2vg.o
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

test : hal2vg
	cd tests && prove -v small.t
