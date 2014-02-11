SRC_COMMON = ../src_common/
SRC_MPI = ../src_mpi/
BOINC_DIR = ../../..
BOOST_DIR = ../..
BOINC_API_DIR = $(BOINC_DIR)/api
BOINC_LIB_DIR = $(BOINC_DIR)/lib
BOINC_ZIP_DIR = $(BOINC_DIR)/zip
FREETYPE_DIR = /usr/include/freetype2

CPP = g++
CPPFLAGS = -O2 \
    -Wall -W -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fno-common \
    -DAPP_GRAPHICS -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS -D NDEBUG \
    -I$(BOINC_DIR) \
	-I$(BOOST_DIR) \
    -I$(BOINC_LIB_DIR) \
    -I$(BOINC_API_DIR) \
    -I$(BOINC_ZIP_DIR) \
    -I$(FREETYPE_DIR) \
    -I$(SRC_COMMON) \
    -I$(SRC_MPI) \
    -L /usr/X11R6/lib \
    -L.

release: pdsat_crypto

libstdc++.a:
	ln -s `g++ -print-file-name=libstdc++.a`

ttfont.cpp:
	ln -s ../../api/ttfont.cpp .

clean: distclean

distclean:
	/bin/rm -f $(PROGS) *.o libstdc++.a

pdsat_crypto: minisat22_wrapper.o mpi_base.o addit_func.o System.o Options.o Solver.o pdsat_crypto.o libstdc++.a $(BOINC_API_DIR)/libboinc_api.a $(BOINC_LIB_DIR)/libboinc.a
	$(CPP) $(CPPFLAGS) minisat22_wrapper.o mpi_base.o addit_func.o System.o Options.o Solver.o pdsat_crypto.o -o pdsat_crypto libstdc++.a -pthread \
	$(BOINC_API_DIR)/libboinc_api.a \
	$(BOINC_LIB_DIR)/libboinc.a

pdsat_crypto.o: ./pdsat_crypto.cpp
	${CPP} ${CPPFLAGS} ./pdsat_crypto.cpp -c

Options.o: ../src_common/minisat/utils/Options.cc
	${CPP} ${CPPFLAGS} ../src_common//minisat/utils/Options.cc -c

System.o: ../src_common/minisat/utils/System.cc
	${CPP} ${CPPFLAGS} ../src_common/minisat/utils/System.cc -c

Solver.o: ../src_common//minisat/core/Solver.cc
	${CPP} ${CPPFLAGS} ../src_common//minisat/core/Solver.cc -c

minisat22_wrapper.o: ../src_common/minisat22_wrapper.cpp
	${CPP} ${CPPFLAGS} ../src_common/minisat22_wrapper.cpp -c

mpi_base.o: ../src_mpi/mpi_base.cpp
	${CPP} ${CPPFLAGS} ../src_mpi/mpi_base.cpp -c
	
addit_func.o: ../src_common/addit_func.cpp
	${CPP} ${CPPFLAGS} ../src_common/addit_func.cpp -c