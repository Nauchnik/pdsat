
#ifndef dminisat_h
#define dminisat_h

#ifdef _WIN32
#define inline __inline // compatible with MS VS
#endif

#include "dminisat_solver.h"

#define UINT_LEN 32
#define FULL_MASK_LEN 26
#define SOLVING_TIME_LEN 15 // info about time of solving tasks
#define MAX_ASSIGNS_COUNT 800

//=================================================================================================
// Public interface:

extern lbool_dminisat parse_DIMACS( FILE * in, solver* s );

extern int dminisat_solve( char* CNFfile, unsigned int full_mask[FULL_MASK_LEN], 
						   unsigned int part_mask[FULL_MASK_LEN], unsigned int values[FULL_MASK_LEN], 
						   int solver_type, int core_len, double corevars_activ_type, int* process_sat_count, 
						   int** b_SAT_set_array, int sort_type, double* cnf_time_from_node,
						   int IsPredict, int IsHardProblem );

#endif
