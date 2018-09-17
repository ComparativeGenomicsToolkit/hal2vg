# simplest possible to start.  dangerous since all header deps done manually. ne
rootPath = ./
include ./include.mk

sidegraphInc = ${sgExportPath}/sidegraph.h ${sgExportPath}/sgcommon.h ${sgExportPath}/sgsequence.h ${sgExportPath}/sgposition.h ${sgExportPath}/sgside.h ${sgExportPath}/sgjoin.h ${sgExportPath}/sgsegment.h ${sgExportPath}/sglookup.h ${hal2sgPath}/sgbuilder.h ${sgExportPath}/side2seq.h

all : hal2vg

cleanFast : 
	rm -f hal2vg hal2vg.o sg2vgproto.o

clean :
	rm -f hal2vg hal2vg.o
	cd deps/sonLib && make clean
	cd deps/hal && make clean
	cd deps/hal2sg && make clean

sg2vgproto.o : sg2vgproto.cpp sg2vgproto.h ${sidegraphInc} ${basicLibsDependencies}
	${cpp} ${cppflags} -I . sg2vgproto.cpp -c

hal2vg.o : hal2vg.cpp sg2vgproto.h ${sidegraphInc} ${basicLibsDependencies}
	${cpp} ${cppflags} -I . hal2vg.cpp -c

${sonLibPath}/sonLib.a :
	cd deps/sonLib && make

${halPath}/halLib.a : ${sonLibPath}/sonLib.a
	cd deps/hal && make

${hal2sgPath}/libhal2sg.a : ${halPath}/halLib.a
	cd deps/hal2sg && make

hal2vg : hal2vg.o sg2vgproto.o ${basicLibsDependencies}
	cd deps/hal2sg && make
	${cpp} ${cppflags} -pthread hal2vg.o sg2vgproto.o ${basicLibs}  -o hal2vg

test : hal2vg
	cd tests && VGDIR=${PWD}/${VGDIR} prove -v small.t
