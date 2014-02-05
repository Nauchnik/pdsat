CPP = g++
SRC_COMMON = ../src_common/
BOOST = ../../
MYSQL = /usr/include/mysql
INCDIRS  = -I${SRC_COMMON} -I${MYSQL} -I${BOOST}

CPP_DEBUG_FLAGS = -O0 -DDEBUG_ALL -g ${INCDIRS}
CPP_RELEASE_FLAGS = -O2 -D NDEBUG -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS ${INCDIRS} 
DCAPI_MASTER_FLAGS = -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -I/usr/include/uuid
DCAPI_MASTER_LIBS = -ldcapi-boinc -lglib-2.0 -luuid /usr/lib/libmysqlclient.so

#all: release debug

release: CPPFLAGS = $(CPP_RELEASE_FLAGS)
release: pdsat-master

debug: CPPFLAGS = $(CPP_DEBUG_FLAGS)
debug: pdsat-master

pdsat-master: master.o minisat22_wrapper.o latin_squares.o addit_func.o Solver.o System.o Options.o
	${CPP} ${CPPFLAGS} master.o minisat22_wrapper.o latin_squares.o addit_func.o Solver.o System.o Options.o -o pdsat-master ${DCAPI_MASTER_LIBS} ${DCAPI_MASTER_FLAGS}
	
Options.o: ../src_common/minisat/utils/Options.cc
	${CPP} ${CPPFLAGS} ../src_common/minisat/utils/Options.cc -c
	
System.o: ../src_common/minisat/utils/System.cc
	${CPP} ${CPPFLAGS} ../src_common/minisat/utils/System.cc -c

Solver.o: ../src_common/minisat/core/Solver.cc
	${CPP} ${CPPFLAGS} ../src_common/minisat/core/Solver.cc -c

minisat22_wrapper.o: ../src_common/minisat22_wrapper.cpp
	${CPP} ${CPPFLAGS} ../src_common/minisat22_wrapper.cpp -c

latin_squares.o: ../src_common/latin_squares.cpp 
	${CPP} ${CPPFLAGS} ../src_common/latin_squares.cpp -c

addit_func.o: ../src_common/addit_func.cpp 
	${CPP} ${CPPFLAGS} ../src_common/addit_func.cpp -c

master.o: master.cpp
	${CPP} ${CPPFLAGS} master.cpp -c

clean:
	rm -fr *.o
	rm pdsat-master*
	clear
	