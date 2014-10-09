CPP = mpiicpc
INC1 = ../../boost/
INC2 = ../src_common/
INC3 = ../src_common/minisat/
INCDIRS = -I${INC1} -I${INC2} -I${INC3}
CPPFLAGS = -O3 -std=c++0x ${INCDIRS} -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS -D _MPI

pd-sat: pd-sat.o Bivium.o mpi_base.o mpi_solver.o mpi_predicter.o addit_func.o minisat22_wrapper.o Solver.o System.o Options.o
	${CPP} ${CPPFLAGS} pd-sat.o Bivium.o mpi_base.o mpi_solver.o mpi_predicter.o addit_func.o minisat22_wrapper.o Solver.o System.o Options.o -o pd-sat
	
Options.o: ../src_common/minisat/utils/Options.cc
	${CPP} ${CPPFLAGS} ../src_common/minisat/utils/Options.cc -c
	
System.o: ../src_common/minisat/utils/System.cc
	${CPP} ${CPPFLAGS} ../src_common/minisat/utils/System.cc -c

Solver.o: ../src_common/minisat/core/Solver.cc
	${CPP} ${CPPFLAGS} ../src_common/minisat/core/Solver.cc -c
	
minisat22_wrapper.o: ../src_common/minisat22_wrapper.cpp
	${CPP} ${CPPFLAGS} ../src_common/minisat22_wrapper.cpp -c

addit_func.o: ../src_common/addit_func.cpp
	${CPP} ${CPPFLAGS} ../src_common/addit_func.cpp -c

mpi_predicter.o: mpi_predicter.cpp 
	${CPP} ${CPPFLAGS} mpi_predicter.cpp -c	

mpi_solver.o: mpi_solver.cpp 
	${CPP} ${CPPFLAGS} mpi_solver.cpp -c
	
mpi_base.o: mpi_base.cpp 
	${CPP} ${CPPFLAGS} mpi_base.cpp -c

Bivium.o: ../src_common/Bivium.cpp 
	${CPP} ${CPPFLAGS} ../src_common/Bivium.cpp -c
	
pd-sat.o: pd-sat.cpp
	${CPP} ${CPPFLAGS} pd-sat.cpp -c
	
clean:
	rm -rf *.o
	rm pd-sat
	clear