CC = mpiicc
CPP = mpiicpc
INC1 = ../../boost/
INC2 = ../src_common/
INC3 = ../src_common/minisat2-2_Jan2013/
INCDIRS = -I${INC1} -I${INC2} -I${INC3}
CFLAGS = -O3 -std=c99 -lz -lm -Wall -Wmissing-prototypes -Wstrict-prototypes -Wsign-compare -D NDEBUG -D _MPI
CPPFLAGS = -O3 ${INCDIRS} -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS -D _MPI

pd-sat_conseq: main.o latin_squares.o addit_func.o minisat22_wrapper.o Solver.o System.o Options.o
	${CPP} ${CPPFLAGS} main.o latin_squares.o addit_func.o minisat22_wrapper.o Solver.o System.o Options.o -o pd-sat_conseq

Options.o: ../src_common/minisat2-2_Jan2013/minisat/utils/Options.cc
	${CPP} ${CPPFLAGS} ../src_common/minisat2-2_Jan2013/minisat/utils/Options.cc -c
	
System.o: ../src_common/minisat2-2_Jan2013/minisat/utils/System.cc
	${CPP} ${CPPFLAGS} ../src_common/minisat2-2_Jan2013//minisat/utils/System.cc -c
	
Solver.o: ../src_common/minisat2-2_Jan2013/minisat/core/Solver.cc
	${CPP} ${CPPFLAGS} ../src_common/minisat2-2_Jan2013/minisat/core/Solver.cc -c
	
minisat22_wrapper.o: ../src_common/minisat22_wrapper.cpp
	${CPP} ${CPPFLAGS} ../src_common/minisat22_wrapper.cpp -c

addit_func.o: ../src_common/addit_func.cpp
	${CPP} ${CPPFLAGS} ../src_common/addit_func.cpp -c

latin_squares.o: ../src_common/latin_squares.cpp
	${CPP} ${CPPFLAGS} ../src_common/latin_squares.cpp -c
	
main.o: main.cpp
	${CPP} ${CPPFLAGS} main.cpp -c
	
clean:
	rm -rf *.o
	rm pd-sat
	clear