CPP = g++
INC1 = /usr/include/mysql
INC2 = /var/lib/boinc/pdsat/src/pdsat/src_common/
INCDIRS = -I${INC1} -I${INC2}
CPPLIBS = /usr/lib/libmysqlclient.so
CPPFLAGS = -O2 ${INCDIRS}

boinc_pdsat_assimilator: boinc_pdsat_assimilator.o Mols.o
	${CPP} ${CPPFLAGS} ${CPPLIBS} boinc_pdsat_assimilator.o Mols.o -o boinc_pdsat_assimilator

boinc_pdsat_assimilator.o: boinc_pdsat_assimilator.cpp
	${CPP} ${CPPFLAGS} boinc_pdsat_assimilator.cpp -c

Mols.o: /var/lib/boinc/pdsat/src/pdsat/src_common/Mols.cpp
	${CPP} ${CPPFLAGS} /var/lib/boinc/pdsat/src/pdsat/src_common/Mols.cpp -c

clean:
	rm -rf *.o
	rm boinc_pdsat_assimilator
	clear