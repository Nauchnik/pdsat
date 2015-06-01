CPP = g++
INCDIRS = -I/usr/include/mysql
CPPLIBS = /usr/lib/libmysqlclient.so
CPPFLAGS = -O2 ${INCDIRS}

boinc_pdsat_assimilator: boinc_pdsat_assimilator.o
	${CPP} ${CPPFLAGS} ${CPPLIBS} boinc_pdsat_assimilator.o -o boinc_pdsat_assimilator

boinc_pdsat_assimilator.o: boinc_pdsat_assimilator.cpp
	${CPP} ${CPPFLAGS} boinc_pdsat_assimilator.cpp -c

clean:
	rm -rf *.o
	rm boinc_pdsat_assimilator
	clear