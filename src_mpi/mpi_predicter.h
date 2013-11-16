// class for predicting of SAT-problem solve time using MPI
#ifndef mpi_predicter_h
#define mpi_predicter_h

#include "mpi_base.h"
#include "../src_common/addit_func.h"

const int	 MAX_STEP_CNF_IN_SET		   = 30;

const int    MAX_WORD_LENGTH			   = 64;
const int    MAX_LINE_LENGTH               = 524288;
const int    MAX_LINE_LENGTH_2             = 8192;
const int    MEDIUM_STRING_LEN             = 256;
const double TRANSP_COAST                  = 0.000001;
const int    NUM_KEY_BITS                  = 64;
const int    MAX_VAR_FOR_RANDOM            = 60;
const double MIN_STOP_TIME				   = 0.01;
const int    MAX_DISTANCE_TO_RECORD        = 20;

struct unchecked_area
{
	boost::dynamic_bitset<> center; // point - center of area. i here means variable # i
	boost::dynamic_bitset<> checked_points; // i here corresponds to checked point with component # i
	int radius;
	int center_count; // count of 1s in center
};

struct checked_area
{
	boost::dynamic_bitset<> center; // point - center of area
	int radius;
	int center_count;
};

struct decomp_set
{
	vector<int> var_choose_order; // indexes of variables in decomp set 
	int set_var_count; // count of known vars - may be different for every set
	int cur_var_changing; // how many vars were changed in current decomp_set relatively to last best point
	bool IsAddedToL2;
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

	// deep predict params
	int deep_diff_decomp_set_count;
	int max_var_deep_predict;
	int deep_predict;
	// how many new best points were finded with such count of new vars
	vector<int> global_count_var_changing; // how many points were found with particular Hamming distance
	int deep_predict_cur_var;
	bool IsRestartNeeded;
	bool IsDecDecomp;
	bool IsSimulatedGranted;
	bool IsFirstPoint;
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
	string deep_predict_file_name;
	string var_activity_file_name;
	vector< vector<int> > combinations;
	list<checked_area> L1; // areas where all points were checked
	list<unchecked_area> L2; // areas where not all points were checked
	unchecked_area current_unchecked_area; // for creating list of points for checking
	// TODO ? vector of unchecked areas where vector of checked must be changed cause of last point 
	int ts_strategy;
	unsigned max_sat_problems;
	double current_predict_start_time;
	double current_predict_time;
	double whole_deep_time;
	double whole_get_deep_tasks_time;
	double whole_get_predict_time;
	double whole_add_new_unchecked_area_time;
	int predict_every_sec;
	double start_sample_varianse_limit;
	
	int cnf_in_set_count;
	int decomp_set_count;
	int total_decomp_set_count;
	int solved_tasks_count;
	string predict_file_name;
	int record_count;
	bool IsFirstStage;
	bool IsRecordUpdated;
	unsigned max_L2_hamming_distance;

	double *var_activity;
	vector<double> total_var_activity;

	int slow_cnf_mask_index;
	vector<decomp_set> decomp_set_arr;
	vector< vector<unsigned> > part_mask_arr;
	vector< vector<unsigned> > all_values_arr;
	vector<double> cnf_real_time_arr;
	vector<char> cnf_status_arr;

	// array of block sum lengths for mass predict
	// in fact it is count of vars for paralleling
	vector<int> sum_block_lens_arr;
	vector<int> sorted_index_array;

	bool MPI_Predict( int argc, char **argv );
	bool ControlProcessPredict( int ProcessListNumber, stringstream &sstream_control );
	bool ComputeProcessPredict( );
	bool GetPredict( );

	bool DeepPredictMain( );
	bool DeepPredictFindNewUncheckedArea( stringstream &sstream );
	bool GetDeepPredictTasks( );
	void GetInitPoint( );
	void NewRecordPoint( int set_index );
	void DecVarLinerDeepMode( stringstream &sstream, fstream &deep_predict_file, int &ProcessListNumber );
	void GetDecDeepPredictTasks( );
	bool IsPointInCheckedArea( boost::dynamic_bitset<> &point );
	bool IsPointInUnCheckedArea( boost::dynamic_bitset<> &point );
	void AddNewUncheckedArea( boost::dynamic_bitset<> &point, stringstream &sstream );
	
	void AllocatePredictArrays( int &cur_tasks_count );
	void ChangeVarChooseOrder( vector<int> var_choose_order, unsigned cur_vars_changing, 
		                       unsigned current_var_count, vector<int> &new_var_choose_order );
	void GetNewHammingPoint( vector<int> var_choose_order, int change_var_count, int &current_var_count, 
		                     vector<int> diff_vec, vector<int> &new_var_choose_order );

	bool PrepareForPredict( );
	bool GetStandartPredictTasks( );

	void GetValuesBasedOnPartMask( unsigned shortcnf_count, int &all_values_arr_index );
	bool GetRandomValuesArray( unsigned shortcnf_count, vector< vector<unsigned> > &values_arr );
	bool IfSimulatedGranted( double predict_time );
	bool WritePredictToFile( int all_skip_count, double whole_time_sec );

private:
	int best_var_num;
	double best_predict_time;

	int real_best_var_num;
	double real_best_predict_time;
	double best_predict_time_last_area; // best time from previous area oc points
	vector<int> real_var_choose_order;

	vector<unsigned> set_len_arr; // array of indexes of sets
	vector<int> cnf_to_stop_arr;  
	vector<int> set_index_arr;
	vector<char> set_status_arr;
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