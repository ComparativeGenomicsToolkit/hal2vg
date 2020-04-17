# simplest possible to start.  dangerous since all header deps done manually. ne
rootPath = ./
include ./include.mk

sidegraphInc = ${sgExportPath}/sidegraph.h ${sgExportPath}/sgcommon.h ${sgExportPath}/sgsequence.h ${sgExportPath}/sgposition.h ${sgExportPath}/sgside.h ${sgExportPath}/sgjoin.h ${sgExportPath}/sgsegment.h ${sgExportPath}/sglookup.h ${hal2sgPath}/sgbuilder.h ${sgExportPath}/side2seq.h

all : hal2vg

cleanFast : 
	rm -f hal2vg hal2vg.o sg2vghandle.o

clean :
	rm -f hal2vg hal2vg.o
	cd deps/sonLib && make clean
	cd deps/hal && make clean
	cd deps/hal2sg && make clean
	cd deps/libbdsg-easy && make clean

sg2vghandle.o : sg2vghandle.cpp sg2vghandle.h ${sidegraphInc} ${basicLibsDependencies}
	${cpp} ${cppflags} -I . sg2vghandle.cpp -c

hal2vg.o : hal2vg.cpp sg2vghandle.h ${sidegraphInc} ${basicLibsDependencies}
	${cpp} ${cppflags} -I . hal2vg.cpp -c

${sonLibPath}/sonLib.a :
	cd deps/sonLib && make

${halPath}/halLib.a : ${sonLibPath}/sonLib.a
	cd deps/hal && make

${hal2sgPath}/libhal2sg.a : ${halPath}/halLib.a
	cd deps/hal2sg && make

${libbdsgPath}/libbdsg.a :
	cd deps/libbdsg-easy && make

${libbdsgPath}/lib/libhandlegraph.a : ${libbdsgPath}/libbdsg.a

${libbdsgPath}/lib/libsdsl.a : ${libbdsgPath}/libbdsg.a

${libbdsgPath}/lib/libdivsufsort.a : ${libbdsgPath}/libbdsg.a

${libbdsgPath}/lib/libdivsufsort64.a : ${libbdsgPath}/libbdsg.a

hal2vg : hal2vg.o sg2vghandle.o ${basicLibsDependencies}
	cd deps/hal2sg && make
	${cpp} ${cppflags} -fopenmp -pthread hal2vg.o sg2vghandle.o ${basicLibs}  -o hal2vg

test : hal2vg
	cd tests && prove -v small.t
