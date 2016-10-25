# simplest possible to start.  dangerous since all header deps done manually. ne
rootPath = ./
include ./include.mk

sidegraphInc = ${sgExportPath}/sidegraph.h ${sgExportPath}/sgcommon.h ${sgExportPath}/sgsequence.h ${sgExportPath}/sgposition.h ${sgExportPath}/sgside.h ${sgExportPath}/sgjoin.h ${sgExportPath}/sgsegment.h ${sgExportPath}/sglookup.h ${hal2sgPath}/sgbuilder.h ${sg2vgPath}/sg2vgjson.h

all : hal2vg

clean : 
	rm -f hal2vg hal2vg.o

cleanAll :
	rm -f hal2vg hal2vg.o
	cd deps/sonLib && make clean
	cd deps/hal && make clean
	cd deps/hal2sg && make clean
	cd deps/sg2vg && make clean

hal2vg.o : hal2vg.cpp  ${sidegraphInc} ${basicLibsDependencies}
	${cpp} ${cppflags} -I . hal2vg.cpp -c

hal2vg :  hal2vg.o ${basicLibsDependencies}
	${cpp} ${cppflags} hal2vg.o ${basicLibs}  -o hal2vg 

