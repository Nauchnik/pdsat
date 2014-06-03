
#ifndef dminisat_h
#define dminisat_h

#ifdef _WIN32
#define inline __inline // compatible with MS VS
#endif

#include "dminisat_solver.h"

#define MAX_ASSIGNS_COUNT 128

//=================================================================================================
// Public interface:

extern lbool_dminisat parse_DIMACS( FILE * in, solver* s );

#endif
