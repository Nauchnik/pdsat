// class for predicting of SAT-problem solve time using MPI
#ifndef mpi_predicter_h
#define mpi_predicter_h

#include "mpi_base.h"
#include "Bivium.h"
#include <thread>

#ifndef WIN32
#include <unistd.h>.
#include <sys/types.h>
#endif

const unsigned MAX_STEP_CNF_IN_SET		   = 30;
const int      MAX_WORD_LENGTH			   = 64;
const int      MAX_LINE_LENGTH             = 524288;
const int      MAX_LINE_LENGTH_2           = 8192;
const int      MEDIUM_STRING_LEN           = 256;
const int      NUM_KEY_BITS                = 64;
const int      MAX_VAR_FOR_RANDOM          = 60;
const double   MIN_STOP_TIME			   = 0.01;
const int      MAX_DISTANCE_TO_RECORD      = 20;
const int      PREDICT_TIMES_COUNT         = 6;
const int      TS2_POINTS_COUNT            = 100;
const unsigned MAX_POW_VALUE               = 1000;
const double   MIN_PERCENT_SOLVED_IN_TIME  = 1.0;

struct unchecked_area
{
	boost::dynamic_bitset<> center; // point - center of area. i here means variable # i
	boost::dynamic_bitset<> checked_points; // i here corresponds to checked point with component # i
	int radius;
	double med_var_activity;
	bool is_partly_checked; // if checked in window mode, then not all points in radius are checked, but we can't chosse it again
};

struct checked_area
{
	boost::dynamic_bitset<> center; // point - center of area
	int radius;
};

struct decomp_set
{
	std::vector<int> var_choose_order; // indexes of variables in decomp set 
	//int set_var_count; // count of known vars - may be different for every set
	int cur_var_changing; // how many vars were changed in current decomp_set relatively to last best point
	bool IsAddedToL2;
	double med_var_activity;
	double diff_variable_activity; // variable by which set differs from current center point
};

struct point_struct
{
    unsigned ones_count;
    unsigned size;
};

class MPI_Predicter : public MPI_Base
{
public:
	// Constructor/Destructor:
    MPI_Predicter( );
    ~MPI_Predicter( );

	int predict_from;
	int predict_to;
	int proc_count;
	int block_count;
	std::string evaluation_type;

	// deep predict params
	int deep_diff_decomp_set_count;
	int max_var_deep_predict;
	int deep_predict;
	// how many new best points were finded with such count of new vars
	std::vector<int> global_count_var_changing; // how many points were found with particular Hamming distance
	int deep_predict_cur_var;
	bool IsRestartNeeded;
	bool IsDecDecomp;
	bool isSimulatedGranted;
	bool isFirstPoint;
	int stop_message;

	double cur_temperature;
	double min_temperature;
	double temperature_multiply_koef;
	double start_temperature_koef; // multiply this koef to predict time and agin temperature
	double point_admission_koef; // how exactly worse can new point be. all others will be interrupted
	double delta;
	double exp_value;
	int cur_vars_changing;
	int global_deep_point_index;
	int global_checked_points_count;
	int global_stopped_points_count;
	int global_skipped_points_count;
	int decomp_sets_in_block;
	std::string deep_predict_file_name;
	std::string var_activity_file_name;
	std::vector< std::vector<int> > combinations;
	std::list<checked_area> L1; // areas where all points were checked
	std::list<unchecked_area> L2; // areas where not all points were checked
	unchecked_area current_unchecked_area; // for creating list of points for checking
	// TODO ? vector of unchecked areas where vector of checked must be changed cause of last point 
	unsigned ts_strategy;
	double current_predict_start_time;
	double current_predict_time;
	double whole_deep_time;
	double whole_get_deep_tasks_time;
	double whole_get_predict_time;
	double whole_add_new_unchecked_area_time;
	int predict_every_sec;
	unsigned unupdated_count;
	double start_sample_variance_limit;
	double prev_area_best_predict_time;
	double predict_time_limit_step;
	
	Problem cnf;
	//unsigned prev_best_decomp_set_power;
	//unsigned prev_best_sum;
	unsigned blob_var_count; // max count of var in decompositions set for writing blob
	unsigned cur_point_number;
	int isSolvedOnPreprocessing;
	std::string tmp_cnf_process_name;
	std::string current_cnf_out_name;
	long long template_cnf_size;
	std::stringstream template_sstream;
	bool IsFirstStage;
	unsigned max_L2_hamming_distance;
	
	std::vector< std::vector<bool> > stream_vec_vec;
	std::vector< std::vector<bool> > state_vec_vec;
	std::vector< std::string > oneliteral_string_vec;
	
	bool MPI_Predict( int argc, char **argv );
	bool ControlProcessPredict( int ProcessListNumber, std::stringstream &sstream_control );
	bool ComputeProcessPredict();
	bool GetPredict();
	bool solverProgramCalling( Minisat::vec<Minisat::Lit> &dummy );
	bool solverSystemCalling( Minisat::vec<Minisat::Lit> &dummy );
	double getCurPredictTime( unsigned cur_var_num, int cur_cnf_in_set_count, unsigned i );
	
	bool DeepPredictMain( );
	bool DeepPredictFindNewUncheckedArea( std::stringstream &sstream );
	bool GetDeepPredictTasks();
	void GetInitPoint();
	void NewRecordPoint( int set_index );
	bool IsPointInCheckedArea( boost::dynamic_bitset<> &point );
	bool IsPointInUnCheckedArea( boost::dynamic_bitset<> &point );
	void AddNewUncheckedArea( boost::dynamic_bitset<> &point, std::stringstream &sstream );
	
	void AllocatePredictArrays();
	
	bool PrepareForPredict();
	bool GetRandomValuesArray( unsigned shortcnf_count,std:: vector< std::vector<unsigned> > &values_arr );
	bool checkSimulatedGranted( double predict_time );
	bool WritePredictToFile( int all_skip_count, double whole_time_sec );
	void SendPredictTask( int ProcessListNumber, unsigned process_number_to_send, int &cur_task_index, 
		                  unsigned &cur_decomp_set_index );
	std::vector<int> BitsetToIntVecPredict( boost::dynamic_bitset<> &bs );
	boost::dynamic_bitset<> IntVecToBitsetPredict( std::vector<int> &variables_vec );
private:
	std::vector<decomp_set> decomp_set_arr;
	std::vector<double> cnf_real_time_arr;
	std::vector<int> cnf_not_solved_check_count;
	std::vector<bool> cnf_issat_arr;
	std::vector<int> cnf_prepr_arr;
	std::vector<int> cnf_status_arr;
	std::vector<double> total_var_activity;
	// array of block sum lengths for mass predict. in fact it is count of vars for paralleling
	std::vector<int> sum_block_lens_arr;
	std::vector<int> sorted_index_array;
	std::vector<bool> vec_IsDecompSetSendedToProcess;
	std::vector<double> predict_time_limites;
	
	int best_var_num;
	long double best_predict_time;
	long double best_predict_time_arr[PREDICT_TIMES_COUNT];
	long double best_sum_time;
	int best_cnf_in_set_count;
	unsigned best_solved_in_time;
	long double best_time_limit;
	
	unsigned total_decomp_set_count;
	unsigned solved_tasks_count;
	std::string predict_file_name;
	unsigned record_count;
	bool IsRecordUpdated;
	double *var_activity;

	int *array_message;
	unsigned array_message_size;

	int real_best_var_num;
	long double real_best_predict_time;
	long double best_predict_time_last_area; // best time from previous area of points
	std::vector<int> real_var_choose_order;
	
	std::vector<unsigned> set_len_arr; // array of indexes of sets
	std::vector<unsigned> cnf_to_stop_arr;  
	std::vector<unsigned > set_index_arr;
	std::vector<int> set_status_arr;
	std::vector<unsigned> node_list;
	std::vector<int> stopped_cnf_count_arr;				
	std::vector<int> skipped_cnf_count_arr; 
	std::vector<int> solved_cnf_count_arr;	
	std::vector<double> cnf_start_time_arr;		    
	std::vector<double> cnf_final_time_arr; 
	std::vector<double> sum_time_arr;
	std::vector<unsigned> solved_in_time_arr;
	std::vector<double> time_limit_arr; // time limit with best predict value for a point in RoEsTe mode
	std::vector<long double> med_time_arr;
	std::vector<long double> predict_time_arr;
	std::vector<long double> predict_part_time_arr;
};

#endif