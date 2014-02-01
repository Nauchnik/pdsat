SRC_COMMON = ../src_common/
MINISAT = ../src_common/minisat2-2_Jan2013/
BOINC_DIR = ../../..
BOINC_API_DIR = $(BOINC_DIR)/api
BOINC_LIB_DIR = $(BOINC_DIR)/lib
BOINC_ZIP_DIR = $(BOINC_DIR)/zip
FREETYPE_DIR = /usr/include/freetype2

CPP = g++
CPPFLAGS = -O3 \
    -Wall -W -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fno-common \
    -DAPP_GRAPHICS -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS -D NDEBUG \
    -I$(BOINC_DIR) \
    -I$(BOINC_LIB_DIR) \
    -I$(BOINC_API_DIR) \
    -I$(BOINC_ZIP_DIR) \
    -I$(FREETYPE_DIR) \
    -I$(SRC_COMMON) \
    -I$(MINISAT) \
    -L /usr/X11R6/lib \
    -L.

release: pdsat

libstdc++.a:
	ln -s `g++ -print-file-name=libstdc++.a`

ttfont.cpp:
	ln -s ../../api/ttfont.cpp .

clean: distclean

distclean:
	/bin/rm -f $(PROGS) *.o libstdc++.a

pdsat: minisat22_wrapper.o latin_squares.o addit_func.o System.o Options.o Solver.o pdsat-client.o libstdc++.a $(BOINC_API_DIR)/libboinc_api.a $(BOINC_LIB_DIR)/libboinc.a
	$(CPP) $(CPPFLAGS) minisat22_wrapper.o latin_squares.o addit_func.o System.o Options.o Solver.o pdsat-client.o -o pdsat libstdc++.a -pthread \
	$(BOINC_API_DIR)/libboinc_api.a \
	$(BOINC_LIB_DIR)/libboinc.a

pdsat-client.o: ./pdsat-client.cpp
	${CPP} ${CPPFLAGS} ./pdsat-client.cpp -c

Options.o: ../src_common/minisat2-2_Jan2013/minisat/utils/Options.cc
	${CPP} ${CPPFLAGS} ../src_common/minisat2-2_Jan2013/minisat/utils/Options.cc -c

System.o: ../src_common/minisat2-2_Jan2013/minisat/utils/System.cc
	${CPP} ${CPPFLAGS} ../src_common/minisat2-2_Jan2013/minisat/utils/System.cc -c

Solver.o: ../src_common/minisat2-2_Jan2013/minisat/core/Solver.cc
	${CPP} ${CPPFLAGS} ../src_common/minisat2-2_Jan2013/minisat/core/Solver.cc -c

minisat22_wrapper.o: ../src_common/minisat22_wrapper.cpp
	${CPP} ${CPPFLAGS} ../src_common/minisat22_wrapper.cpp -c

latin_squares.o: ../src_common/latin_squares.cpp
	${CPP} ${CPPFLAGS} ../src_common/latin_squares.cpp -c

addit_func.o: ../src_common/addit_func.cpp
	${CPP} ${CPPFLAGS} ../src_common/addit_func.cpp -c