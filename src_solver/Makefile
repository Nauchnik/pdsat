CC = gcc
CPP = g++
INC1 = ../src_common/
INC2 = ../src_common/minisat22/
INCDIRS  = -I${INC1} -I${INC2}
# -D NDEBUG to skip assert(), -W... are warnings options
CCFLAGS = -std=c99 -lz -Wall -Wmissing-prototypes -Wstrict-prototypes -Wsign-compare
#CPPFLAGS =  
#CC_COMMON_FLAGS= -std=c99 -lz -Wall -Wmissing-prototypes -Wstrict-prototypes -Wsign-compare
CC_DEBUG_FLAGS= -O0 -D DEBUG_ALL -g 
CC_RELEASE_FLAGS = -O2 -D NDEBUG
#CPP_COMMON_FLAGS = -O2
CPP_DEBUG_FLAGS = -O0 ${INCDIRS} -DDEBUG_ALL -g
CPP_RELEASE_FLAGS = -O2 ${INCDIRS}
LIBS_DEBUG = -pthread  -lcrypto -lstdc++ -lm -Wl
LIBS_RELEASE = -static -pthread -lstdc++ -lm -Wl

#all: release debug

#release: clean
release: CCFLAGS  += $(CC_RELEASE_FLAGS)
release: CPPFLAGS = $(CPP_RELEASE_FLAGS)
release: minisat22_mod

#debug: clean
debug: CPPFLAGS = $(CPP_DEBUG_FLAGS)
#debug: pdsat-master_DEBUG pdsat-client_DEBUG
debug: minisat22_mod_DEBUG

minisat22_mod: main.o minisat22_wrapper.o latin_squares.o Solver.o System.o Options.o
	${CPP} ${CPPFLAGS} main.o minisat22_wrapper.o latin_squares.o Solver.o System.o Options.o -o minisat22_mod ${LIBS_RELEASE}
	
Options.o: ../src_common/minisat22/utils/Options.cc
	${CPP} ${CPPFLAGS} ../src_common/minisat22/utils/Options.cc -c
	
System.o: ../src_common/minisat22/utils/System.cc
	${CPP} ${CPPFLAGS} ../src_common/minisat22/utils/System.cc -c
	
Solver.o: ../src_common/minisat22/core/Solver.cc
	${CPP} ${CPPFLAGS} ../src_common/minisat22/core/Solver.cc -c

latin_squares.o: ../src_common/latin_squares.cpp 
	${CPP} ${CPPFLAGS} ../src_common/latin_squares.cpp -c
	
minisat22_wrapper.o: ../src_common/minisat22/minisat22_wrapper.cpp
	${CPP} ${CPPFLAGS} ../src_common/minisat22/minisat22_wrapper.cpp -c
	          
client.o: main.cpp
	${CPP} ${CPPFLAGS} main.cpp -c

clean:
 clean_obj
	rm minisat22_mod*
	clear

clean_obj:
	rm -fr *.o