// class for predicting of SAT-problem solve time using MPI
#ifndef mpi_predicter_h
#define mpi_predicter_h

#include "mpi_base.h"
#include "Bivium.h"

const unsigned MAX_STEP_CNF_IN_SET		   = 30;

const int    MAX_WORD_LENGTH			   = 64;
const int    MAX_LINE_LENGTH               = 524288;
const int    MAX_LINE_LENGTH_2             = 8192;
const int    MEDIUM_STRING_LEN             = 256;
const double TRANSP_COAST                  = 0.000001;
const int    NUM_KEY_BITS                  = 64;
const int    MAX_VAR_FOR_RANDOM            = 60;
const double MIN_STOP_TIME				   = 0.01;
const int    MAX_DISTANCE_TO_RECORD        = 10;
const int    PREDICT_TIMES_COUNT            = 6;

struct unchecked_area
{
	boost::dynamic_bitset<> center; // point - center of area. i here means variable # i
	boost::dynamic_bitset<> checked_points; // i here corresponds to checked point with component # i
	int radius;
	double med_var_activity;
	vector<double> predict_times;
	double cur_er_predict_time; // for current er - for comparison in sorting
	double sum_time; // sum of predict time 
};

struct checked_area
{
	boost::dynamic_bitset<> center; // point - center of area
	int radius;
};

struct decomp_set
{
	vector<int> var_choose_order; // indexes of variables in decomp set 
	//int set_var_count; // count of known vars - may be different for every set
	int cur_var_changing; // how many vars were changed in current decomp_set relatively to last best point
	bool IsAddedToL2;
	double med_var_activity;
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
	string evaluation_type;

	// deep predict params
	int deep_diff_decomp_set_count;
	int max_var_deep_predict;
	int deep_predict;
	int er_strategy;
	// how many new best points were finded with such count of new vars
	vector<int> global_count_var_changing; // how many points were found with particular Hamming distance
	int deep_predict_cur_var;
	bool IsRestartNeeded;
	bool IsDecDecomp;
	bool IsSimulatedGranted;
	bool isFirstPoint;

	double cur_temperature;
	double min_temperature;
	double temperature_multiply_koef;
	double start_temperature_koef; // multiply this koef to predict time and agin temperature
	double point_admission_koef; // how exactly worse can new point be. all others will be interrupted
	double delta;
	double exp_value;
	double exp_denom;
	int cur_vars_changing;
	int global_deep_point_index;
	int global_checked_points_count;
	int global_stopped_points_count;
	int global_skipped_points_count;
	int decomp_sets_in_block;
	string deep_predict_file_name;
	string var_activity_file_name;
	vector< vector<int> > combinations;
	list<checked_area> L1; // areas where all points were checked
	list<unchecked_area> L2; // areas where not all points were checked
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
	//unsigned prev_best_decomp_set_power;
	//unsigned prev_best_sum;
	unsigned max_var_count_state_writing;

	vector< vector<bool> > stream_vec_vec;
	vector< vector<bool> > state_vec_vec;
	
	bool MPI_Predict( int argc, char **argv );
	bool ControlProcessPredict( int ProcessListNumber, stringstream &sstream_control );
	bool ComputeProcessPredict();
	bool GetPredict();
	
	bool DeepPredictMain( );
	bool DeepPredictFindNewUncheckedArea( stringstream &sstream );
	bool GetDeepPredictTasks();
	void GetInitPoint();
	void NewRecordPoint( int set_index );
	bool IsPointInCheckedArea( boost::dynamic_bitset<> &point );
	bool IsPointInUnCheckedArea( boost::dynamic_bitset<> &point );
	void AddNewUncheckedArea( boost::dynamic_bitset<> &point, vector<double> &cur_predict_times, double sum_time, stringstream &sstream );
	
	void AllocatePredictArrays();
	
	bool PrepareForPredict();
	bool GetRandomValuesArray( unsigned shortcnf_count, vector< vector<unsigned> > &values_arr );
	bool IfSimulatedGranted( double predict_time );
	bool WritePredictToFile( int all_skip_count, double whole_time_sec );
	void SendPredictTask( int ProcessListNumber, int process_number_to_send, int &cur_task_index, unsigned &cur_decomp_set_index );
	vector<int> BitsetToIntVecPredict( boost::dynamic_bitset<> &bs );
	boost::dynamic_bitset<> IntVecToBitsetPredict( vector<int> &variables_vec );
	
	bool IsFirstStage;
	unsigned max_L2_hamming_distance;

private:
	vector<decomp_set> decomp_set_arr;
	vector<double> cnf_real_time_arr;
	vector<bool> cnf_issat_arr;
	vector<int> cnf_prepr_arr;
	vector<int> cnf_status_arr;
	vector<double> total_var_activity;
	// array of block sum lengths for mass predict. in fact it is count of vars for paralleling
	vector<int> sum_block_lens_arr;
	vector<int> sorted_index_array;
	vector<bool> vec_IsDecompSetSendedToProcess;
	
	int best_var_num;
	double best_predict_time;
	double best_predict_time_arr[PREDICT_TIMES_COUNT];
	double best_sum_time;
	int best_cnf_in_set_count;
	int *array_message; // for sending via MPI
	unsigned array_message_size;
	
	unsigned total_decomp_set_count;
	unsigned solved_tasks_count;
	string predict_file_name;
	unsigned record_count;
	bool IsRecordUpdated;
	double *var_activity;

	int real_best_var_num;
	double real_best_predict_time;
	double best_predict_time_last_area; // best time from previous area of points
	vector<int> real_var_choose_order;
	
	vector<unsigned> set_len_arr; // array of indexes of sets
	vector<int> cnf_to_stop_arr;  
	vector<unsigned > set_index_arr;
	vector<int> set_status_arr;
	vector<int> node_list;
	vector<int> stopped_cnf_count_arr;				
	vector<int> skipped_cnf_count_arr; 
	vector<int> solved_cnf_count_arr;	
	vector<double> cnf_start_time_arr;		    
	vector<double> cnf_final_time_arr; 
	vector<double> sum_time_arr;			    
	vector<double> med_time_arr;
	vector<double> predict_time_arr;
	vector<double> predict_part_time_arr;
};

#endif