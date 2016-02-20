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

const int      MAX_CORE_LEN = 10000;
const double   MIN_SOLVE_TIME = 0.0000001;
const unsigned RECOMMEND_BATCH_VAR_COUNT = 21;
const unsigned MAX_BATCH_VAR_COUNT = 31;
const unsigned MAX_PART_MASK_VAR_COUNT = 24;
const unsigned UINT_LEN = 32;
const unsigned FULL_MASK_LEN = 120;
const unsigned SOLVING_TIME_LEN = 15; // info about time of solving tasks
const double   SOLVER_PARSE_SIMP_TIME = 0.03; // solver parse time + simplification time
const long double HUGE_DOUBLE = 1e+308;

enum ProblemStates{Solved, SolvedOnPreprocessing, Interrupted};

class MPI_Base
{
protected:
	std::vector< std::vector<int> > clause_array;
	std::vector<int> clause_lengths;
	std::vector<bool> b_SAT_set_array;
	std::vector< std::vector<unsigned> > values_arr;
	std::map<int, unsigned> core_var_indexes; // indexes of variables in core set
	boost::random::mt19937 gen;
	unsigned known_vars_count;
public:
    MPI_Base();
    ~MPI_Base();
	
	int rank;
	int corecount;
	std::string solver_name;
	unsigned koef_val;
	std::string schema_type;
	unsigned core_len;
	double start_activity;
	bool isConseq;
	bool isPredict;
	bool isMakeSatSampleAnyWay;
	
	unsigned full_mask[FULL_MASK_LEN];
	unsigned part_mask[FULL_MASK_LEN];
	unsigned mask_value[FULL_MASK_LEN]; // particular value of bits which set in part_mask
	
	unsigned activity_vec_len;
	std::string known_point_file_name;
	bool isSolverSystemCalling; // calling of solver file by system command

	// Common CNF input data
	unsigned var_count;
	unsigned lit_count;
	unsigned clause_count;
	unsigned full_mask_var_count;
	unsigned part_mask_var_count;
	unsigned all_tasks_count;
	double te;
	unsigned first_stream_var_index;
	unsigned known_last_bits;
	unsigned keystream_len;
	unsigned cnf_in_set_count;
	unsigned input_var_num;
	int current_task_index; 
	int process_sat_count;
	double cnf_time_from_node;
	bool isPlainText;
	std::string evaluation_type;

	int sat_count;
	int verbosity;
	bool isSolveAll;
	double max_solving_time; // max time in seconds for solving particular SAT subproblem
	int max_nof_restarts;
	std::string input_cnf_name;
	std::vector<int> var_choose_order;
	std::vector<int> full_var_choose_order; // all variables that can be chosen to decomp set

	// Read header "p cnf [var_count] [clause_count]" from DIMACS file
	bool ReadVarCount( );
	// Read clauses from DIMACS file
	bool ReadIntCNF( );
	
	// Make array var_choose_order with vars sorted by given rule
	bool MakeVarChoose( );
	bool makeStandardMasks( unsigned &part_var_power );
	bool makeIntegerMasks(std::vector<std::vector<int>> cartesian_elements);
	bool GetMainMasksFromVarChoose( std::vector<int> &var_choose_order );
	bool GetValuesFromVarChoose( unsigned &part_var_power );
	bool getValuesFromIntegers(std::vector<std::vector<int>> cartesian_elements);
	void MakeSatSample( std::vector< std::vector<bool> > &state_vec_vec, 
		                std::vector< std::vector<bool> > &stream_vec_vec,
						std::vector< std::vector<bool> > &addit_vec_vec,
						int rank = 0 );
	void MakeSingleSatSample(
			std::vector<bool> &state_vec, 
			std::vector<bool> &stream_vec, 
			const int seed_num,
			Solver *S,
			std::vector<lbool> predefined_vars = std::vector<lbool>());
	bool AnalyzeSATset( double cnf_time_from_node );
	bool CheckSATset( std::vector<int> &lit_SAT_set_array );
	bool MakeAssignsFromMasks( unsigned *full_mask, 
							   unsigned *part_mask, 
						       unsigned *value, 
							   Minisat::vec< Minisat::vec<Minisat::Lit> > &dummy_vec );
	//bool MakeAssignsFromFile( int current_task_index, unsigned long long before_binary_length, Minisat::vec< Minisat::vec<Minisat::Lit> > &dummy_vec );
	
	void MakeRandArr( std::vector< std::vector<unsigned> > &rand_arr, unsigned shortcnf_count, unsigned rnd_uint32_count );
	void MakeUniqueRandArr( std::vector<unsigned> &rand_arr, unsigned rand_arr_len, unsigned max_rand_val );
	std::string MakeSolverLaunchString( std::string solver_name, std::string cnf_name, double maxtime_seconds_str );
};

#endif
