CC = gcc
CPP = g++
# -D NDEBUG to skip assert(), -W... are warnings options
CCFLAGS = -std=c99 -lz -Wall -Wmissing-prototypes -Wstrict-prototypes -Wsign-compare
#CPPFLAGS =  
#CC_COMMON_FLAGS= -std=c99 -lz -Wall -Wmissing-prototypes -Wstrict-prototypes -Wsign-compare
CC_DEBUG_FLAGS= -O0 -D DEBUG_ALL -g 
CC_RELEASE_FLAGS = -O2 -D NDEBUG
#CPP_COMMON_FLAGS = -O2
CPP_DEBUG_FLAGS = -O0 -DDEBUG_ALL -g
CPP_RELEASE_FLAGS = -O2 
DCAPI_CLIENT_LIBS_DEBUG = -pthread -ldcapi-client-boinc -lboinc_api -lboinc -ldcapi-boinc -lcrypto -lstdc++ -lm -Wl
DCAPI_CLIENT_LIBS_RELEASE = -static -pthread -ldcapi-client-boinc -lboinc_api -lboinc -ldcapi-boinc -lcrypto -lstdc++ -lm -Wl
DCAPI_MASTER_CFLAGS = -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -I/usr/include/uuid
DCAPI_MASTER_LIBS = -ldcapi-boinc -lglib-2.0 -luuid

#all: release debug

#release: clean
release: CCFLAGS  += $(CC_RELEASE_FLAGS)
release: CPPFLAGS = $(CPP_RELEASE_FLAGS)
release: assimilator

#debug: clean
debug: CCFLAGS  += $(CC_DEBUG_FLAGS)
debug: CPPFLAGS = $(CPP_DEBUG_FLAGS)
debug: pdsat-master_DEBUG pdsat-client_DEBUG

assimilator: master.o mpi_base.o
	${CPP} ${CPPFLAGS} master.o mpi_base.o -o assimilator ${DCAPI_MASTER_LIBS} ${DCAPI_MASTER_CFLAGS}

pdsat-client: client.o dminisat.o dminisat_solver.o mpi_base.o
	${CPP} ${CPPFLAGS} client.o dminisat.o dminisat_solver.o mpi_base.o -o pdsat-client ${DCAPI_CLIENT_LIBS_RELEASE}

pdsat-master_DEBUG: master.o mpi_base.o
	${CPP} ${CPPFLAGS} master.o mpi_base.o -o pdsat-master_DEBUG ${DCAPI_MASTER_LIBS} ${DCAPI_MASTER_CFLAGS}

pdsat-client_DEBUG: client.o dminisat.o dminisat_solver.o mpi_base.o
	${CPP} ${CPPFLAGS} client.o dminisat.o dminisat_solver.o mpi_base.o -o pdsat-client_DEBUG ${DCAPI_CLIENT_LIBS_DEBUG}

dminisat.o: dminisat.c
	${CPP} ${CPPFLAGS}  dminisat.c -c
	
dminisat_solver.o: dminisat_solver.c
	${CPP} ${CCFLAGS} dminisat_solver.c -c

mpi_base.o: mpi_base.cpp 
	${CPP} ${CPPFLAGS} mpi_base.cpp -c		
	
master.o: master.cpp
	${CPP} ${CPPFLAGS} master.cpp -c
	
client.o: client.cpp
	${CPP} ${CPPFLAGS} client.cpp -c

clean: clean_obj
	rm assimilator*
	clear

clean_obj:
	rm -fr *.o