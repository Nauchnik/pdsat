// class for solving SAT-problem using MPI

#ifndef mpi_base_h
#define mpi_base_h

#include <cstdio>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>

/*#ifdef __cplusplus
extern "C" 
{
#endif
#include "dminisat.h"
#ifdef __cplusplus
}
#endif*/

#include "dminisat.h"

#ifdef _WIN32
#include "..\include_win_pdsat\dirent.h"
#else
#include <dirent.h>
#endif

//#include "minisat2.h"

#define UINT_LEN 32
#define FULL_MASK_LEN 34
#define MAX_CORE_LEN 1024

using namespace std;

class MPI_Base
{
public:
	// Constructor/Destructor:
    MPI_Base( );
    ~MPI_Base( );

	int corecount;
	int solver_type;    
	int koef_val;
	int schema_type; 
	int core_len;
	double corevars_activ_type;
	bool IsConseq;
	int check_every_conflict;
	
	unsigned int full_mask[FULL_MASK_LEN];
	unsigned int part_mask[FULL_MASK_LEN];

	// PB input data
	int constr_clauses_count;
	int obj_clauses_count;
	int obj_vars_count;
	int obj_vars[40]; // boolean view of int value

	// Common CNF input data
	int var_count;
	int lit_count;
	int clause_count;

	int full_mask_var_count;
	int part_mask_var_count;

	int all_tasks_count;

	bool IsSAT;
	bool IsPB; //  pseudo Boolean mode. if 0 then common CNF mode
	int PB_mode; // Pseudo Boolean mode. 1 - inequality mode, 2 - equality mode
	int best_lower_bound;
	int upper_bound;
	int verbosity;

	char *input_cnf_name;

	int *b_SAT_set_array;
	int **clause_array;
	int *clause_lengths;
	int **lits_clause_array;
	int *lits_clause_lengths;
	int *lit_SAT_set_array;
	
	int var_choose_order[MAX_CORE_LEN];
	
	// Read header "p cnf [var_count] [clause_count]" from DIMACS file
	bool ReadVarCount( );
	// Read clauses from DIMACS file
	bool ReadIntCNF( );
	
	// Make array var_choose_order with vars sorted by given rule
	bool MakeVarChoose( );
	bool MakeStandartMasks( unsigned long long int &part_var_power, 
							unsigned int **&values_arr );
	bool GetMainMasksFromVarChoose( );
	bool GetValuesFromVarChoose( unsigned long long int &part_var_power, 
								 unsigned int **&values_arr );
	
	bool AnalyzeSATset( );
	bool cpuTimeInHours( double full_seconds, int &real_hours, int &real_minutes, 
		                 int &real_seconds );
	bool WriteTimeToFile( double whole_time_sec );
	bool CheckSATset( int &int_answer );
	bool SortValuesDecrease( unsigned int range, unsigned int *&sorted_index_array );

	// service procedures
	void shl64( unsigned long long int &val_for_left_shift, unsigned int bit_count );
	void equalize_arr( unsigned int arr1[FULL_MASK_LEN], unsigned int arr2[FULL_MASK_LEN] );
	int getdir( string dir, vector<string> &files );
};

#endif
