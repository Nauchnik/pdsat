CPP = g++
INCDIRS = -I/usr/include/mysql
CPPLIBS = /usr/lib/libmysqlclient.so
CPPFLAGS = -O2 ${INCDIRS}

boinc_assimilator: boinc_assimilator.o
	${CPP} ${CPPFLAGS} ${CPPLIBS} boinc_assimilator.o -o boinc_assimilator

boinc_assimilator.o: boinc_assimilator.cpp
	${CPP} ${CPPFLAGS} boinc_assimilator.cpp -c

clean_obj:
	rm -fr *.o