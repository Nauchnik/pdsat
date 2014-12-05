CPP = g++
BOOST = ../../
MYSQL = /usr/include/mysql
UUID = /usr/include/uuid
INCDIRS  = -I${MYSQL} -I${BOOST} -I${UUID}

CPP_DEBUG_FLAGS = -O0 -DDEBUG_ALL -g ${INCDIRS}
CPP_RELEASE_FLAGS = -O2 -D NDEBUG -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS ${INCDIRS} -MYSQL
MYSQL_LIBS = -luuid /usr/lib/libmysqlclient.so

#all: release debug

release: CPPFLAGS = $(CPP_RELEASE_FLAGS)
release: pdsat_work_generator

debug: CPPFLAGS = $(CPP_DEBUG_FLAGS)
debug: pdsat_work_generator

pdsat_work_generator: work_generator.o
	${CPP} ${CPPFLAGS} work_generator.o -o pdsat_work_generator $(MYSQL_LIBS)

work_generator.o: work_generator.cpp
	${CPP} ${CPPFLAGS} work_generator.cpp -c
	
clean:
	rm -fr *.o
	rm pdsat_work_generator*
	clear
	