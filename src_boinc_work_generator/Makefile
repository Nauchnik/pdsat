CPP = g++
BOOST = ../../
MYSQL = /usr/include/mysql
COMMON = ../src_common/
INCDIRS  = -I${MYSQL} -I${BOOST} -I${COMMON}

CPP_DEBUG_FLAGS = -O0 -DDEBUG_ALL -g ${INCDIRS}
CPP_RELEASE_FLAGS = -O2 -D NDEBUG -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS ${INCDIRS} -MYSQL
DCAPI_MASTER_FLAGS = -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -I/usr/include/uuid
DCAPI_MASTER_LIBS = -ldcapi-boinc -lglib-2.0 -luuid /usr/lib/libmysqlclient.so

#all: release debug

release: CPPFLAGS = $(CPP_RELEASE_FLAGS)
release: pdsat-crypto-master

debug: CPPFLAGS = $(CPP_DEBUG_FLAGS)
debug: pdsat-crypto-master

pdsat-crypto-master: master.o addit_func.o
	${CPP} ${CPPFLAGS} master.o addit_func.o -o pdsat-crypto-master ${DCAPI_MASTER_LIBS} ${DCAPI_MASTER_FLAGS}

master.o: master.cpp
	${CPP} ${CPPFLAGS} master.cpp -c

addit_func.o: ../src_common/addit_func.cpp
	${CPP} ${CPPFLAGS} ../src_common/addit_func.cpp -c
	
clean:
	rm -fr *.o
	rm pdsat-crypto-master*
	clear
	