// +---------------------------------------------------------------------------+
// | PDSAT: Parallel and Distributed SAT solver                                |
// | MPI-realization of packet D-SAT for SAT-problem solving                   |                         
// +---------------------------------------------------------------------------+
// | Institute for System Dynamics and Control Theory SB RAS                   |   
// +---------------------------------------------------------------------------+
// | Author: Oleg Zaikin <zaikin.icc@gmail.com>                                |
// +---------------------------------------------------------------------------+

#include "mpi_base.h"
#include "mpi_solver.h"
#include "mpi_predicter.h"

struct Flags
{
	unsigned predict_to;
	std::string solver_name;
	int koef_val;
	std::string schema_type; 
	int proc_count; 
	int core_len; 
	bool IsPredict;	
	int cnf_in_set_count;
	double start_activity; 			
	bool IsConseq;				
	int verbosity;
	int deep_predict;
	int max_var_deep_predict;
	double start_temperature_koef;
	double point_admission_koef;
	int ts_strategy;
	bool IsFirstStage;
	int max_L2_hamming_distance;
	bool IsSolveAll;
	double max_solving_time;
	int max_nof_restarts;
	int skip_tasks;
	std::string evaluation_type;
	double max_solving_time_koef;
	bool no_increm;
	double te; // for (ro es te) predict strategy
	int variables_each_integer;
};

// prototypes
//---------------------------------------------------------
bool GetInputFlags( int &argc, char **&argv, Flags &myflags );
const char* hasPrefix( const char* str, const char* prefix );
std::string hasPrefix_String( std::string str, std::string prefix );
void WriteUsage( );
void TestPredict( );
void TestDeepPredict( );
void TestSolve( );

//---------------------------------------------------------
int main( int argc, char** argv )
{
// main procedure
	int full_mask_var_count;
	char *input_cnf_name;
	
#ifdef _DEBUG
	TestSolve();
	TestDeepPredict();
	//TestSolveQAP();
	//TestPBSolve();
	//TestSATSolve();
#endif

	Flags myflags;

	// check input flags
	if ( !GetInputFlags( argc, argv, myflags ) ) { 
		printf( "\n Error in GetInputFlags" ); 
		WriteUsage( ); 
		return 1; 
	}
	
	if ( myflags.IsConseq ) {
		MPI_Solver mpi_s;
		mpi_s.solver_name = myflags.solver_name;
		if ( myflags.core_len != -1 )
			mpi_s.core_len = myflags.core_len;
		mpi_s.verbosity = myflags.verbosity;
		mpi_s.MPI_ConseqSolve( argc, argv );
		std::cout << std::endl << "Correct end of MPI_ConseqSolve";
		return 0;
	}

	if ( !(myflags.IsConseq) ) {
		// check count of input parameters
		if ( argc < 2 ) {
			std::cout << "\n argc must be >= 2" << std::endl;
			WriteUsage( );
			return 1;
		}
		input_cnf_name = argv[1]; // get name of file with input CNF
	}
	
	if ( myflags.IsPredict ) {
		MPI_Predicter mpi_p;
		mpi_p.input_cnf_name          = input_cnf_name; // for C dminisat
		mpi_p.solver_name             = myflags.solver_name;
		if ( myflags.koef_val != -1 )
			mpi_p.koef_val            = myflags.koef_val;
		if ( myflags.schema_type != "" )
			mpi_p.schema_type         = myflags.schema_type;
		if ( myflags.proc_count != -1 )
			mpi_p.proc_count          = myflags.proc_count;
		if ( myflags.core_len != -1 )
			mpi_p.core_len            = myflags.core_len;
		if ( myflags.cnf_in_set_count != -1 )
			mpi_p.cnf_in_set_count    = myflags.cnf_in_set_count;
		if ( myflags.start_activity != -1 )
			mpi_p.start_activity = myflags.start_activity;
		if ( myflags.verbosity > 0 )
			mpi_p.verbosity = myflags.verbosity;
		if ( myflags.max_var_deep_predict > 0 )
			mpi_p.max_var_deep_predict = myflags.max_var_deep_predict;
		if ( myflags.deep_predict > 0 )
			mpi_p.deep_predict = myflags.deep_predict;
		if ( myflags.start_temperature_koef > 0 )
			mpi_p.start_temperature_koef = myflags.start_temperature_koef;
		if ( myflags.point_admission_koef > 0 )
			mpi_p.point_admission_koef = myflags.point_admission_koef;
		if ( myflags.ts_strategy != -1 )
			mpi_p.ts_strategy = myflags.ts_strategy;
		if ( myflags.max_L2_hamming_distance > 0 )
			mpi_p.max_L2_hamming_distance = myflags.max_L2_hamming_distance;
		if ( myflags.max_solving_time > 0 )
			mpi_p.max_solving_time = myflags.max_solving_time;
		if ( myflags.evaluation_type != "" )
			mpi_p.evaluation_type = myflags.evaluation_type;
		if ( myflags.predict_to > 0 )
			mpi_p.predict_to = myflags.predict_to;
		mpi_p.IsFirstStage = myflags.IsFirstStage;
		if ( myflags.te > 0 )
			mpi_p.te = myflags.te;

		// Predict compute cost
		if ( !mpi_p.MPI_Predict( argc, argv ) ) {
			printf( "\n Error in MPI_Predict" );
			return 1;
		}
	}
	else { // Solve SAT-problem
		// get count of vars for splitting
		if ( argc == 3 ) {
			full_mask_var_count = atoi( argv[2] );
			if ( full_mask_var_count < 0 ) {
				printf( "\n Error. full_mask_var_count < 0" );
				WriteUsage( );
				return 1;
			}
		}
		else 
			full_mask_var_count = 31; // default

		MPI_Solver mpi_s;
		mpi_s.input_cnf_name = input_cnf_name;
		
		if ( full_mask_var_count != -1 )
			mpi_s.full_mask_var_count = (unsigned)full_mask_var_count;

		if ( myflags.solver_name != "" )
			mpi_s.solver_name = myflags.solver_name;
		if ( myflags.koef_val != -1 )
			mpi_s.koef_val    = myflags.koef_val;
		if ( myflags.schema_type != "" )
			mpi_s.schema_type = myflags.schema_type;
		if ( myflags.core_len != -1 )
			mpi_s.core_len    = myflags.core_len;
		if ( myflags.start_activity != -1 )
			mpi_s.start_activity = myflags.start_activity;
		if ( myflags.max_solving_time > 0 )
			mpi_s.max_solving_time = myflags.max_solving_time;
		if ( myflags.max_nof_restarts > 0 )
			mpi_s.max_nof_restarts = myflags.max_nof_restarts;
		if ( myflags.max_solving_time_koef > 0 )
			mpi_s.max_solving_time_koef = myflags.max_solving_time_koef;
		if (myflags.evaluation_type != "")
			mpi_s.evaluation_type = myflags.evaluation_type;
		if (myflags.te > 0)
			mpi_s.te = myflags.te;
		if (myflags.variables_each_integer > 0 )
			mpi_s.variables_each_integer = myflags.variables_each_integer;
		mpi_s.no_increm = myflags.no_increm;

		mpi_s.verbosity		   = myflags.verbosity;
		mpi_s.isConseq         = myflags.IsConseq;
		mpi_s.isSolveAll       = myflags.IsSolveAll;

		if ( !mpi_s.MPI_Solve( argc, argv ) ) {
			printf( "\n Error in MPI_Solve" );
			return 1;
		}
		mpi_s.~MPI_Solver( );
	}

	return 0;
}

//---------------------------------------------------------
void WriteUsage( )
{
// Write info about usage
	std :: cout << 
	"\n USAGE: pdsat [options] <input_cnf> <split_var_count>"
	"\n options::"
	"\n   -solver = { minisat, <path_to_solver_file> }, -s = -solver"
	"\n		minisat - minisat 2.2"
	"\n   -koef = { 1, 2, ... }"
	"\n		coefficient for tasks count calculating"
	"\n	  -schema = { first, missone, firstandlast }"
	"\n		0 = first        - 1, 2, 3, ..."
	"\n		1 = missone      - 1, 3, 5, ..."
	"\n		2 = literal count"
	"\n		3 = jeroslaw-wang"
	"\n     4 = implicant count"
	"\n     a5 = 31 vars 1..9, 20..30, 42..52"
	"\n     rslos_end = last cells of of rslos"
	"\n   -cnf_in_set    - count of CNFs in set"
	"\n   -proc_count    - count of cores for predicting"
	"\n   -core          - count of core vars for dm and m2mod (by default = 64)"
	"\n   -conseq        - conseq solve of cnf in current path"
	"\n   -file_assumptions = read assumptions from file"
	"\n   -verb = {0, 1, 2, 3}  - scale of output info"
	"\n   -solve_all     - solve all problems, even if SAT was found"
    "\n   *** deep predict ***"
	"\n   -deep_predict - (integer) mode of deep predict. 1 - fixed neighborhood of sets for every current best point. "
	"\n			2 - new neighborhood for every new best point"
	"\n   -max_var_deep - (float) koefficient (0,1). koef*predict_to = max amount of variables for changing"
    "\n   -ts_strategy - way of choosing new unchecked area";
	"\n   -max_sat_problems - max count of sat problems in iteration of prediction";
	"\n   -max_L2_hamming_distance - max distance of point from L2 from current point";
	"\n   -no_first_stage - disable first stage mode (checking all areas until power of dec. set increase)";
	"\n   -max_solving_time - max time in seconds for splving particular SAT subproblem";
	"\n   -max_nof_restarts - maximum number of restarts";
	"\n   -skip_tasks - count of already solved tasks";
	"\n   -eval - type of prediction evaluation {time, propagation}";
	"\n   -no_increm - disable incremental mode while solving";
	"\n   -te - time limit for backdoor strategy in predict";
}

//---------------------------------------------------------
bool hasPrefix_String( std::string str, std::string prefix, std::string &value )
{
	int found = str.find( prefix );
	if ( found != -1 ) {
		value = str.substr( found + prefix.length( ) );
		return true;
	}
	return false;
}

//---------------------------------------------------------
bool GetInputFlags( int &argc, char **&argv, Flags &myflags )
{
// Get input keys if such exists 
	int i, k;
	bool IsSchemaFixed = false;
	std::stringstream sstream;
	std::string argv_string, value;
	// default values
	myflags.solver_name		    = "";
	myflags.koef_val			= -1,
	myflags.schema_type			= "";
	myflags.proc_count			= -1;
	myflags.start_activity      = -1;
	myflags.core_len			= -1;
	myflags.IsConseq            = false;
	myflags.verbosity			= -1;
	myflags.core_len            = 0;
	myflags.deep_predict = -1;
	myflags.cnf_in_set_count = -1;
	myflags.start_temperature_koef = -1;
	myflags.point_admission_koef = -1;
	myflags.verbosity = -1;
	myflags.max_var_deep_predict = -1;
	myflags.ts_strategy = -1;
	myflags.IsFirstStage = true;
	myflags.max_L2_hamming_distance = -1;
	myflags.IsSolveAll = false;
	myflags.max_solving_time = 0;
	myflags.max_nof_restarts = 0;
	myflags.skip_tasks = 0;
	myflags.evaluation_type = "";
	myflags.max_solving_time_koef = 0;
	myflags.no_increm = false;
	myflags.te = 0;
	myflags.predict_to = 0;
	myflags.variables_each_integer = 0;
	k = 0;
	
	// check every input parameters for flag existing
    for ( i = 0; i < argc; i++ ) {
		sstream.str( "" );
		sstream.clear();
		sstream << argv[i];
		argv_string = sstream.str( );
		if ( ( hasPrefix_String( argv_string, "-proc=",        value ) ) ||
			      ( hasPrefix_String( argv_string, "-proc_count=",  value ) ) )
		{
			myflags.proc_count = atoi( value.c_str( ) );
			if ( myflags.proc_count < 1 )
			{ std::cerr << "Error. proc_count is negative"; return false; }
		}
		else if ( hasPrefix_String( argv_string, "-conseq", value ) )
			myflags.IsConseq = true;
		else if ( hasPrefix_String( argv_string, "-no_increm", value ) )
			myflags.no_increm = true;
		else if ( hasPrefix_String( argv_string, "-no_first_stage", value ) )
			myflags.IsFirstStage = false;
		else if ( hasPrefix_String( argv_string, "-solve_all", value ) )
			myflags.IsSolveAll = true;
		else if (hasPrefix_String(argv_string, "-integer_variables=", value))
			myflags.variables_each_integer = atoi(value.c_str());
		else if ( hasPrefix_String( argv_string, "-max_L2_hamming_distance=", value ) )
			myflags.max_L2_hamming_distance = atoi( value.c_str( ) );
		else if ( hasPrefix_String( argv_string, "-skip_tasks=", value ) )
			myflags.skip_tasks = atoi( value.c_str() );
		else if (hasPrefix_String(argv_string, "-predict_to=", value))
			myflags.predict_to = atoi(value.c_str());
		else if ( hasPrefix_String( argv_string, "-max_nof_restarts=", value ) )
			myflags.max_nof_restarts = atoi( value.c_str( ) );
		else if ( hasPrefix_String( argv_string, "-max_solving_time=", value ) )
			myflags.max_solving_time = atof( value.c_str( ) );
		else if (hasPrefix_String(argv_string, "-max_solving_time_koef=", value))
			myflags.max_solving_time_koef = atof(value.c_str());
		else if ( hasPrefix_String( argv_string, "-te=", value ) )
			myflags.te = atof( value.c_str( ) );
		//else if ( hasPrefix_String( argv_string, "-penalty=", value ) )
		//	myflags.penalty = atof( value.c_str( ) );
		else if ( hasPrefix_String( argv_string, "-deep_predict=", value ) )  {
			myflags.deep_predict = atoi( value.c_str( ) );
			switch ( myflags.deep_predict ) {
				case 1 : // method 1
					myflags.max_var_deep_predict = 10;
					break;
				case 2 : // method 2
					myflags.max_var_deep_predict = 5;
					break;
				case 3 : // method 3
					myflags.max_var_deep_predict = 3;
					break;
				case 4 : // method 4 : local search with checking all area
					myflags.max_var_deep_predict = 3;
					break;
				case 5 : // method 5 : simulated annelaling
					myflags.max_var_deep_predict = 2; // Hamming distance == 2
					break;
				case 6 : // method 5 : local + simulating
					myflags.max_var_deep_predict = 1; // Hamming distance == 2
					break;
				default :
					std::cout << "***Warning. Incorrect myflags.deep_predict";
					myflags.max_var_deep_predict = 2;
					break;
			}
		}
		else if ( ( hasPrefix_String( argv_string, "-solver=",      value ) ) || 
				  ( hasPrefix_String( argv_string, "-s=",           value ) ) )
		{
			myflags.solver_name = value;
		}
		else if ( hasPrefix_String( argv_string, "-koef=", value ) ) {
			myflags.koef_val = atoi( value.c_str( ) );
			if ( myflags.koef_val < 1 ) {
				myflags.koef_val = 1;
				std::cout << "koef_val was changed to 1";
			}
		}
		else if ( ( hasPrefix_String( argv_string, "-core=",     value ) ) ||
			      ( hasPrefix_String( argv_string, "-core_len=", value ) ) )
		{
			myflags.core_len = atoi( value.c_str( ) );
			if ( myflags.core_len < 0 ) {
				myflags.core_len = 0;
				std::cerr << "Error. core_len < 0, changed to 0";
			}
			if ( myflags.core_len > MAX_CORE_LEN ) { // check < MAX_CORE_LEN
				myflags.core_len = MAX_CORE_LEN;
				std::cerr << "Error. core_len > MAX_CORE_LEN, changed to MAX_CORE_LEN";
			}
		}
		else if ( ( hasPrefix_String( argv_string, "-cnf_in_set=",       value ) ) ||
			      ( hasPrefix_String( argv_string, "-cnf_in_set_count=", value ) ) )
		{
			myflags.cnf_in_set_count = atoi( value.c_str( ) );
			if ( myflags.cnf_in_set_count < 0 ) {
				myflags.cnf_in_set_count = 64;
				std::cout << "cnf_in_set_count was changed to " << myflags.cnf_in_set_count;
			}
			if ( myflags.cnf_in_set_count > ( 1 << MAX_STEP_CNF_IN_SET ) ) {
				myflags.cnf_in_set_count = ( 1 << MAX_STEP_CNF_IN_SET );
				std::cout << "cnf_in_set_count was changed to " << myflags.cnf_in_set_count;
			}
		}
		else if ( hasPrefix_String( argv_string, "-schema=",              value ) ) {
			IsSchemaFixed = true;
			myflags.schema_type = value;
		}
		else if ( hasPrefix_String( argv_string, "-eval=", value ) )
			myflags.evaluation_type = value;
		else if ( hasPrefix_String( argv_string, "-ts_strategy=",          value ) )
			myflags.ts_strategy = atoi( value.c_str( ) );
		else if ( hasPrefix_String( argv_string, "-start_activity=",       value ) )
			myflags.start_activity = atof( value.c_str( ) );
		else if ( hasPrefix_String( argv_string, "-verb=",				   value ) )
			myflags.verbosity = atoi( value.c_str( ) );
		// deep predict
		else if ( hasPrefix_String( argv_string, "-max_var_deep=",		   value ) )
			myflags.max_var_deep_predict = atoi( value.c_str( ) );
		else if ( hasPrefix_String( argv_string, "-start_temperature_koef=", value ) )
			myflags.start_temperature_koef = atof( value.c_str( ) );
		else if ( hasPrefix_String( argv_string, "-point_admission_koef=",	 value ) )
			myflags.point_admission_koef = atof( value.c_str( ) );
		else if ( ( argv_string == "-h"     ) || 
			      ( argv_string == "-help"  ) || 
				  ( argv_string == "--help" ) )
			WriteUsage( );
		else if ( argv_string[0] == '-' )
            std::cerr << "ERROR! unknown flag " << argv_string;
		else
            argv[k++] = argv[i]; // skip flag arguments
    }
    argc = k;
	
	// if all keys needed for predicting exists then start predicting
	if ( ( myflags.deep_predict > 0 ) ||
		 ( myflags.cnf_in_set_count > 0 )
	   )
	{
		myflags.IsPredict = true;
	}
	else myflags.IsPredict = false;

	/*
	cout << "\nproc_count "    << proc_count;
	cout << "\nIsPredict "     << IsPredict;
	*/

	return true;
}

void TestSolve()
{
	Solver S;

	char *input_cnf_name = "../src_common/tresh72_0.cnf";	
	std::string rslos_table_name = "../src_common/bits_4_1/1010.txt";
	int process_sat_count = 0;
	unsigned int full_mask[FULL_MASK_LEN];
	unsigned int part_mask[FULL_MASK_LEN];
	unsigned int value[FULL_MASK_LEN];
	MPI_Solver mpi_s;

	std::vector<std::vector<int>> cartesian_elements;
	cartesian_elements.resize(1);
	cartesian_elements[0].resize(3);
	cartesian_elements[0][0] = 0;
	cartesian_elements[0][1] = 20;
	cartesian_elements[0][2] = 40;
	mpi_s.var_choose_order.resize(60);
	for (unsigned i = 0; i < mpi_s.var_choose_order.size(); i++)
		mpi_s.var_choose_order[i] = i + 1;
	mpi_s.all_tasks_count = cartesian_elements.size();
	mpi_s.values_arr.resize(mpi_s.all_tasks_count);
	for (unsigned i = 0; i < mpi_s.values_arr.size(); ++i) {
		mpi_s.values_arr[i].resize(FULL_MASK_LEN);
		for (unsigned j = 0; j < mpi_s.values_arr[i].size(); ++j)
			mpi_s.values_arr[i][j] = 0;
	}
	mpi_s.makeIntegerMasks(cartesian_elements);

	mpi_s.input_cnf_name = input_cnf_name;
	//mpi_s.schema_type = "rslos_end";
	mpi_s.solver_name = "minisat";
	mpi_s.ReadIntCNF();
	mpi_s.MakeVarChoose();
	int current_task_index = 0;
	
	for ( unsigned i=0; i < FULL_MASK_LEN; i++ )
		full_mask[i] = part_mask[i] = 0;
	full_mask[0] = part_mask[0] = value[0] = 1;
	full_mask[1] = 1048575;
	part_mask[1] = 1023;

	mpi_s.all_tasks_count = 16;
	//mpi_s.verbosity = 2;
	//mpi_s.SolverRun( S, process_sat_count, cnf_time_from_node,current_task_index );
	int core_len = 64;
	double corevars_activ_type = 1;
	int sort_type = 0;
}

//---------------------------------------------------------
void TestDeepPredict( )
{
	std::cout << "*** DEBUG MODE" << std::endl;
	MPI_Predicter mpi_p;

	std::vector<unsigned> rand_arr; 
	unsigned rand_arr_len = 192; 
	unsigned max_rand_val = 56;

	mpi_p.MakeUniqueRandArr( rand_arr, rand_arr_len, max_rand_val );
	sort( rand_arr.begin(), rand_arr.end() );

	mpi_p.core_len = 177;
	mpi_p.cnf_in_set_count = 2000000;
	mpi_p.deep_predict = 6;
	mpi_p.max_var_deep_predict = 1;

	mpi_p.part_mask[0] = 3;
	mpi_p.part_mask[1] = 4846842;
	mpi_p.part_mask[2] = 315856;
	std::vector< std::vector<unsigned> > values_arr;

	mpi_p.GetInitPoint();
	values_arr.resize(10);
	for ( unsigned i = 0; i < 10; i++ )
		values_arr[i].resize(FULL_MASK_LEN);
	mpi_p.GetRandomValuesArray( 10, values_arr );

	mpi_p.schema_type = "bivium_Ending2"; 
	mpi_p.GetInitPoint( );

	unchecked_area ua;
	ua.center.resize( mpi_p.core_len );
	for ( unsigned i = 0; i < mpi_p.core_len; i ++ )
		ua.center.set(i);
	mpi_p.L2.push_back( ua );
	//boost::dynamic_bitset<> bs = IntVecToBitset( mpi_p.core_len, mpi_p.var_choose_order );
	std::stringstream sstream;
	//mpi_p.AddNewUncheckedArea( bs, sstream );
	//mpi_p.current_unchecked_area.center = bs;
	mpi_p.cur_vars_changing = 1;
	mpi_p.isFirstPoint = false;
	//mpi_p.schema_type = "bivium_Ending2"
	mpi_p.GetDeepPredictTasks( );
	mpi_p.cur_vars_changing = 1;
}