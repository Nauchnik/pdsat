// class for solving SAT-problem using MPI

#ifndef mpi_base_h
#define mpi_base_h

#ifdef _MPI
#include <mpi.h>
#endif

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <list>
#include <vector>
#include <map>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <bitset>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/dynamic_bitset.hpp>

#include "minisat22_wrapper.h"
#include "addit_func.h"
using namespace Addit_func;
using namespace std;

const int      MAX_CORE_LEN = 800;
const double   MIN_SOLVE_TIME = 0.000001;
const unsigned RECOMMEND_BATCH_VAR_COUNT = 21;
const unsigned MAX_BATCH_VAR_COUNT = 31;
const unsigned MAX_PART_MASK_VAR_COUNT = 24;
const unsigned UINT_LEN = 32;
const unsigned FULL_MASK_LEN = 26;
const unsigned SOLVING_TIME_LEN = 15; // info about time of solving tasks
const unsigned MAX_ASSIGNS_COUNT = 800;
const double   SOLVER_PARSE_SIMP_TIME = 0.03; // solver parse time + simplification time
const double   HUGE_DOUBLE = 1e+308;

enum ProblemStates{Solved, SolvedOnPreprocessing, Interrupted};

class MPI_Base
{
protected:
	vector< vector<int> > clause_array;
	vector<int> clause_lengths;
	vector<int> b_SAT_set_array;
	vector< vector<unsigned> > values_arr;
	map<int, unsigned> core_var_indexes; // indeces of variables in core set
	boost::random::mt19937 gen;
public:
    MPI_Base();
    ~MPI_Base();
	
	int rank;
	int corecount;
	string solver_name;
	unsigned koef_val;
	string schema_type;
	unsigned core_len;
	double start_activity;
	bool IsConseq;
	int check_every_conflict;
	bool IsPredict;
	bool isMakeSatSampleAnyWay;
	
	unsigned *full_mask;
	unsigned *part_mask;
	unsigned *mask_value; // particular value of bits which set in part_mask
	
	unsigned activity_vec_len;
	string known_point_file_name;
	string base_known_assumptions_file_name;
	string known_assumptions_file_name;
	bool isSolverSystemCalling; // calling of solver file by system command
	
	bool IsPB; //  pseudo Bool mode. if 0 then common CNF mode
	int PB_mode; // Pseudo Boolean mode. 1 - inequality mode, 2 - equality mode
	int best_lower_bound;
	int upper_bound;
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
	unsigned all_tasks_count;
	double te;
	double er;
	double penalty;
	unsigned first_stream_var_index;
	unsigned known_last_bits;
	unsigned keystream_len;
	unsigned cnf_in_set_count;
	unsigned input_var_num;
	int current_task_index; 
	int process_sat_count;
	double cnf_time_from_node;

	int sat_count;
	int verbosity;
	bool IsSolveAll;
	double max_solving_time; // max time in seconds for solving particular SAT subproblem
	unsigned keybit_count;
	unsigned long long assumptions_count;
	string rslos_table_name;
	int max_nof_restarts;
	char *input_cnf_name;
	vector<int> var_choose_order;
	vector<int> full_var_choose_order; // all variables that can be chosen to decomp set
	vector<int> all_vars_set;
	vector<int> rslos_lengths;
	double *all_var_activity;

	// Read header "p cnf [var_count] [clause_count]" from DIMACS file
	bool ReadVarCount( );
	// Read clauses from DIMACS file
	bool ReadIntCNF( );
	
	// Make array var_choose_order with vars sorted by given rule
	bool MakeVarChoose( );
	bool MakeStandardMasks( unsigned &part_var_power );
	bool GetMainMasksFromVarChoose( vector<int> &var_choose_order );
	bool GetValuesFromVarChoose( unsigned &part_var_power );

	void MakeSatSample( vector< vector<bool> > &state_vec_vec, vector< vector<bool> > &stream_vec_vec );
	
	bool AnalyzeSATset( );
	bool CheckSATset( vector<int> &lit_SAT_set_array );
	bool MakeAssignsFromMasks( unsigned *full_mask, 
							   unsigned *part_mask, 
						       unsigned *value, 
							   vec< vec<Lit> > &dummy_vec );
	bool MakeAssignsFromFile( int current_task_index, unsigned long long before_binary_length, vec< vec<Lit> > &dummy_vec );
	
	void MakeRandArr( vector< vector<unsigned> > &rand_arr, unsigned shortcnf_count, unsigned rnd_uint32_count );
	void MakeUniqueRandArr( vector<unsigned> &rand_arr, unsigned rand_arr_len, unsigned max_rand_val );
	std::string make_solver_launch_str( std::string solver_name, std::string cnf_name, double maxtime_seconds_str );
};

#endif
