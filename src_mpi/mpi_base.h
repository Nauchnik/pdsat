// class for solving SAT-problem using MPI

#ifndef mpi_base_h
#define mpi_base_h

#ifdef _MPI
#include <mpi.h>
#endif

#include <cstdio>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <list>
#include <vector>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <bitset>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/dynamic_bitset.hpp>

#ifdef __cplusplus
extern "C" 
{
#endif
#include "dminisat.h"
#ifdef __cplusplus
}
#endif

#ifdef _WIN32
#include "./win_headers/dirent.h"
#else
#include <dirent.h>
#endif

#include "minisat22_wrapper.h"
#include "addit_func.h"
using namespace Addit_func;
using namespace std;

#define MAX_CORE_LEN 800
const double MIN_SOLVE_TIME = 0.000001;
const unsigned MAX_BATCH_VAR_COUNT = 18;
const unsigned MAX_PART_MASK_VAR_COUNT = 24;

enum ProblemStates{ Solved, SolvedOnPreprocessing, Interrupted };

class MPI_Base
{
public:
    MPI_Base( );
    ~MPI_Base( );

	int rank;
	int corecount;
	int solver_type;    
	int koef_val;
	string schema_type; 
	unsigned core_len;
	double start_activity;
	bool IsConseq;
	int check_every_conflict;
	bool IsPredict;
	bool IsFileAssumptions;
	
	unsigned int full_mask[FULL_MASK_LEN];
	unsigned int part_mask[FULL_MASK_LEN];

#ifdef _MPI
	MPI_Datatype mpi_mask;
	MPI_Datatype mpi_solving_time;
	MPI_Datatype mpi_var_activity;
#endif

	unsigned activity_vec_len;
	string known_point_file_name;
	string known_assumptions_file_name;

	// PB input data
	int constr_clauses_count;
	int obj_clauses_count;
	int obj_vars_count;
	int obj_vars[40]; // boolean view of int value

	// Common CNF input data
	unsigned var_count;
	unsigned lit_count;
	unsigned clause_count;

	unsigned full_mask_var_count;
	unsigned part_mask_var_count;

	int all_tasks_count;

	vector<int> rslos_lengths;

	bool sat_count;
	bool IsPB; //  pseudo Boolean mode. if 0 then common CNF mode
	int PB_mode; // Pseudo Boolean mode. 1 - inequality mode, 2 - equality mode
	int best_lower_bound;
	int upper_bound;
	int verbosity;
	int IsHardProblem;
	int sort_type;
	bool IsSolveAll;
	double max_solving_time; // max time in seconds for solving particular SAT subproblem
	unsigned keybit_count;
	unsigned assumptions_string_count;
	string rslos_table_name;
	int max_nof_restarts;
	
	char *input_cnf_name;

	int **clause_array;
	int *clause_lengths;
	int **lits_clause_array;
	unsigned *lits_clause_lengths;
	int *b_SAT_set_array;
	
	vector<int> var_choose_order;
	
	// Read header "p cnf [var_count] [clause_count]" from DIMACS file
	bool ReadVarCount( );
	// Read clauses from DIMACS file
	bool ReadIntCNF( );

	bool SolverRun( Solver *&S, unsigned int *full_mask, unsigned int *part_mask, 
					unsigned int *value, int &process_sat_count, int &current_obj_val,
					double &cnf_time_from_node, double *solving_times,
					int current_task_index );
	
	void AddSolvingTimeToArray( ProblemStates cur_problem_state, double cnf_time_from_node, double *solving_times );
	
	// Make array var_choose_order with vars sorted by given rule
	bool MakeVarChoose( );
	bool MakeStandartMasks( unsigned long long int &part_var_power, 
							unsigned int **&values_arr );
	bool GetMainMasksFromVarChoose( vector<int> &var_choose_order );
	bool GetValuesFromVarChoose( unsigned long long int &part_var_power, 
								 unsigned int **&values_arr );
	
	bool AnalyzeSATset( );
	bool CheckSATset( vector<int> &lit_SAT_set_array );
	void MakeAssignsFromMasks( unsigned full_mask[FULL_MASK_LEN], 
							   unsigned part_mask[FULL_MASK_LEN], 
						       unsigned values[FULL_MASK_LEN], 
							  vec< vec<Lit> > &dummy_vec );
	void MakeAssignsFromFile( int current_task_index, vec< vec<Lit> > &dummy_vec );
	bool SortValuesDecrease( unsigned int range, unsigned int *&sorted_index_array );

	// service procedures
	boost::dynamic_bitset<> IntVecToBitset( unsigned bitset_len, vector<int> &vec_int );
	vector<int> BitsetToIntVec( boost::dynamic_bitset<> &bs );
	void shl64( unsigned long long int &val_for_left_shift, unsigned int bit_count );
	void equalize_arr( unsigned int arr1[FULL_MASK_LEN], unsigned int arr2[FULL_MASK_LEN] );
	void MakeRandArr( vector< vector<unsigned> > &rand_arr, unsigned shortcnf_count, unsigned rnd_uint32_count );
	void MakeUniqueRandArr( vector<unsigned> &rand_arr, unsigned rand_arr_len, 
		              unsigned max_rand_val );
	void PrintVector( vector<int> &vec );

	int getdir( string dir, vector<string> &files );
	unsigned uint_rand();
};

#endif
