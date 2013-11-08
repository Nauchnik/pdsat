
#ifndef dminisat_h
#define dminisat_h

#ifdef _WIN32
#define inline __inline // compatible with MS VS
#endif

/*#ifdef __cplusplus
extern "C"
{
#endif
#include "dminisat_solver.h"
#ifdef __cplusplus
}
#endif*/

#include "dminisat_solver.h"

#define UINT_LEN 32
#define FULL_MASK_LEN 34

//=================================================================================================
// Public interface:

extern int dminisat_solve( char* CNFfile, unsigned int full_mask[FULL_MASK_LEN], 
						   unsigned int part_mask[FULL_MASK_LEN], unsigned int values[FULL_MASK_LEN], 
						   int solver_type, int core_len, double corevars_activ_type, int* IsSATFinded, 
						   int** b_SAT_set_array, int sort_type, double* cnf_time_from_node,
						   int IsPredict, int IsHardProblem, unsigned long long int last_iteration_done,
						   int *max_attempt );

#endif
