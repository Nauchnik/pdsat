#include "mpi_solver.h"

#pragma warning( disable : 4996 )

const int    MAX_WORD_LENGTH			   = 64;
const int    MAX_LINE_LENGTH               = 524288;
const int    MAX_LINE_LENGTH_2             = 8192;
const int    MEDIUM_STRING_LEN             = 256;
const double TRANSP_COAST                  = 0.000001;
const int    NUM_KEY_BITS                  = 64;
//const int MAX_VAR_IN_CNF		        = 65536;
//const int MAX_SIZE_OF_FILE	        = 42949672;
//const int MAX_CNF_IN_FOLDER		    = 65536;
//const int MAX_STEP_CNF_COUNT	        = 30;
//const int LITTLE_STRING_LEN           = 16;

//=============================================================================
// Constructor/Destructor:

MPI_Solver :: MPI_Solver( ) :
	extra_tasks_count      ( 0 ),
	orig_tasks_count       ( 0 ),
	full_mask_tasks_count  ( 0 ),
	exch_activ			   ( 1 ),
	skip_tasks             ( 0 ),
	solving_info_file_name ( "solving_info" )
{ }

MPI_Solver :: ~MPI_Solver( )
{ }

int make_QAP_values( int num_elements, unsigned comb_len, unsigned int **&values_arr, bool IsMakeArray )
{
	int comb_count  = 0;
	int erase_count = 0;
	int real_count  = 0;

	vector<int> elements(num_elements);

	for ( int i = 0; i < num_elements; i++ )
		elements[i] = i;
	
	assert(comb_len > 0 && comb_len <= elements.size());
	vector<unsigned long> positions(comb_len, 0);
	// TODO update with new function for combinations
	/*combinations_recursive( num_elements, comb_len, comb_count, erase_count, 
		                    real_count, elements, comb_len, positions, 0, 0, 
							values_arr, IsMakeArray );*/
	cout << endl << "comb_count is "  << comb_count;
	cout << endl << "erase_count is " << erase_count;
	cout << endl << "real_count is "  << real_count;
	return real_count;
}

//---------------------------------------------------------
bool MPI_Solver :: GetExtraTasks( unsigned int **&values_arr, unsigned int *full_mask_ext )
{
	unsigned int full_mask_temp[FULL_MASK_LEN];

	// index in full_mask_ext that is needed increased by extra 1 vars
	unsigned int extra_addon_index = 0;
	// value that is needed to be added to values with extra 1 vars
	unsigned int extra_addon_value = 0;

	// use full_mask to make it's increased modification
	equalize_arr( full_mask_temp, full_mask );
			
	// temporary action
	full_mask_var_count++;
	part_mask_var_count = full_mask_var_count; // part_mask == full_mask
			
	cout << "\n Start of extra_tasks_count procedures" << endl;
	if ( !MakeVarChoose( ) ) {
		cout << "\n Error in MakeVarChoose" << endl;
		// is't needed to deallocate memory - MPI_Abort will do it
		MPI_Abort( MPI_COMM_WORLD, 0 );
		return 1;
	}
	cout << "\n Correct end of MakeVarChoose" << endl;
			
	if ( !GetMainMasksFromVarChoose( var_choose_order ) ) {
		cout << "\n Error in GetMasksFromVarChoose" << endl;
		return false;
	}
	cout << "\n Correct end of GetMainMasksFromVarChoose" << endl;
	
	// change it back
	full_mask_var_count--;
	part_mask_var_count = full_mask_var_count;
	equalize_arr( full_mask_ext, full_mask );
	equalize_arr( full_mask, full_mask_temp );
	equalize_arr( part_mask, full_mask ); // part_mask == full_mask

	unsigned int shift_count = 0;
	bool IsAddonFound = false;

	extra_addon_index = 1; // init value

	while ( !IsAddonFound && ( extra_addon_index < FULL_MASK_LEN ) )
	{
		// init
		shift_count = 0;
		extra_addon_value = 1;
		extra_addon_index++; // first time it is 2
		while ( !IsAddonFound && ( shift_count < 32 ) )
		{
			if ( ( full_mask[extra_addon_index] & extra_addon_value ) != 
				 ( full_mask_ext[extra_addon_index] & extra_addon_value ) )
			{
				IsAddonFound = true;
				break;
			}
			extra_addon_value <<= 1;
			shift_count++;
		}
	}
	
	cout << "Correct construction of full_mask_ext" << endl;

	if ( !IsAddonFound ) {
		cout << "Error. IsAddonFound == false after GetExtraTasks" << endl;
		return false;
	}

	int i;
	
	// before it last extra_tasks_count values of values_arr are undefined
	// define them basing on first extra_addon_index values
	for( i = 0; i < extra_tasks_count; i++ )
		values_arr[i + orig_tasks_count][extra_addon_index] = 
		values_arr[i + full_mask_tasks_count][extra_addon_index] + extra_addon_value;

	return true;
}

void MPI_Solver :: WriteSolvingTimeInfo( double *solving_times, double *total_solving_times, 
										 unsigned solved_tasks_count, unsigned sat_count, 
										 double finding_first_sat_time )
{
	if ( solving_times[0] < total_solving_times[0] ) // update min time
		total_solving_times[0] = solving_times[0];
	if ( solving_times[1] > total_solving_times[1] ) // update max time
		total_solving_times[1] = solving_times[1];
	total_solving_times[2] += solving_times[2] / all_tasks_count;
	if ( ( total_solving_times[3] == 0 ) && ( solving_times[3] != 0 ) ) // update sat time
		total_solving_times[3] = solving_times[3];

	unsigned long long solved_problems_count = 0;
	for ( unsigned i=4; i < SOLVING_TIME_LEN; i++ ) {
		total_solving_times[i] += solving_times[i];
		solved_problems_count  += (unsigned long long)total_solving_times[i];
	}

	stringstream sstream;
	sstream << "IsFileAssumptions " << IsFileAssumptions << endl;
	sstream << "var_choose_order size " << var_choose_order.size() << endl;
	for ( unsigned i=0; i < var_choose_order.size(); i++ )
		sstream << var_choose_order[i] << " ";

	sstream << endl;
	sstream << "solved_tasks_count     " << solved_tasks_count    << " / " << all_tasks_count << endl;
	sstream << "solved_problems_count  " << solved_problems_count << endl;
	sstream << "*** total_solving_times" << endl;
	sstream << "min                    " << total_solving_times[0] << endl;
	sstream << "max                    " << total_solving_times[1] << endl;
	sstream << "med                    " << total_solving_times[2] << endl;
	sstream << "sat                    " << total_solving_times[3] << endl;
	sstream << "finding_first_sat_time " << finding_first_sat_time << endl;
	sstream << "sat_count              " << sat_count << endl;
	if ( solved_problems_count > 0 ) {
		sstream << "SOLVED BY PREPROCESSING, count ";
		sstream.width( 12 ); 
		sstream << left << fixed << (unsigned)total_solving_times[4]  << " " << total_solving_times[4] * 100 / solved_problems_count << " %" << endl;
		sstream << "(0,   0.0001)   count ";
		sstream.width( 12 ); 
		sstream << left << fixed << (unsigned)total_solving_times[5]  << " " << total_solving_times[5] * 100 / solved_problems_count << " %" << endl;
		sstream << "(0.0001, 0.001) count ";
		sstream.width( 12 ); 
		sstream << left << fixed << (unsigned)total_solving_times[6]  << " " << total_solving_times[6] * 100 / solved_problems_count << " %" << endl;
		sstream << "(0.001, 0.01)   count ";
		sstream.width( 12 );
		sstream << left << fixed << (unsigned)total_solving_times[7]  << " " << total_solving_times[7] * 100 / solved_problems_count << " %" << endl;		
		sstream << "(0.01, 0.1)     count ";
		sstream.width( 12 ); 
		sstream << left << fixed << (unsigned)total_solving_times[8]  << " " << total_solving_times[8] * 100 / solved_problems_count << " %" << endl;
		sstream << "(0.1,    1)     count ";
		sstream.width( 12 ); 
		sstream << left << fixed << (unsigned)total_solving_times[9]  << " " << total_solving_times[9] * 100 / solved_problems_count << " %" << endl;
		sstream << "(1,     10)     count "; 
		sstream.width( 12 ); 
		sstream << left << fixed << (unsigned)total_solving_times[10] << " " << total_solving_times[10] * 100 / solved_problems_count << " %" << endl;
		sstream << "(10,   100)     count "; 
		sstream.width( 12 ); 
		sstream << left << fixed << (unsigned)total_solving_times[11] << " " << total_solving_times[11] * 100 / solved_problems_count << " %" << endl;
		sstream << "(100, 1000)     count ";
		sstream.width( 12 );
		sstream << left << fixed << (unsigned)total_solving_times[12] << " " << total_solving_times[12] * 100 / solved_problems_count << " %" << endl;
		sstream << "(1000, inf)     count "; 
		sstream.width( 12 );
		sstream << left << fixed << (unsigned)total_solving_times[13] << " " << total_solving_times[13] * 100 / solved_problems_count << " %" << endl;
		sstream << "INTERRUPTED, solving time > " << max_solving_time << " ";
		sstream.width( 12 );
		sstream << left << fixed << (unsigned)total_solving_times[14] << " " << total_solving_times[14] * 100 / solved_problems_count << " %" << endl;		
	}
	if ( verbosity > 0 )
		cout << sstream.str();
	ofstream solving_info_file( solving_info_file_name.c_str() );
	solving_info_file << sstream.str();
	solving_info_file.close();
	sstream.clear(); sstream.str("");
}

bool MPI_Solver :: ControlProcessSolve( int first_range_tasks_count, unsigned int *full_mask_ext,
                                        unsigned int **values_arr )
{
	int i;
	unsigned int j;
	int mpi_i = 0;
	int solved_tasks_count = 0;
	int process_sat_count = 0;
	int next_task_index = 0;
	MPI_Status status,
		       current_status;
	unsigned int value[FULL_MASK_LEN];
	unsigned int tmp_mask[FULL_MASK_LEN];
	double start_time = MPI_Wtime(), final_time, finding_first_sat_time = 0;
	
	cout << "\n ControlProcessSolve is running" << endl;

	double total_solving_times[SOLVING_TIME_LEN];
	double solving_times[SOLVING_TIME_LEN];
	total_solving_times[0] = 1 << 30; // start min len
	for ( unsigned i = 1; i < SOLVING_TIME_LEN; i++ )
		total_solving_times[i] = 0;

	for ( i = 0; i < FULL_MASK_LEN; i++ )
		value[i] = 0;

	// send to all cores (except # 0) tasks from 1st range
	for ( int i = skip_tasks; i < skip_tasks + first_range_tasks_count; i++ ) {
		equalize_arr( value, values_arr[i] );
		
		// send new index of task
		MPI_Send( &mpi_i,             1, MPI_INT,  i + 1, 0, MPI_COMM_WORLD );

		// don't send valus when we have file with assimptions
		if ( IsFileAssumptions )
			continue;
		
		MPI_Send( &extra_tasks_count, 1, MPI_INT, i + 1, 0, MPI_COMM_WORLD );

		if ( ( extra_tasks_count ) && ( i >= full_mask_tasks_count ) ) { // modified tasks
			equalize_arr( tmp_mask, full_mask_ext ); // 09.01.2011 checking error
			MPI_Send( &tmp_mask, 1, mpi_mask, mpi_i + 1, 0, MPI_COMM_WORLD );
		}
		else // if no extra tasks or if first orig tasks sending 08.01.11
			MPI_Send( &full_mask,     1, mpi_mask, i + 1, 0, MPI_COMM_WORLD );
	
		// if default mode, part_mask and full_mask can be different
		if ( !extra_tasks_count ) // original tasks
			MPI_Send( &part_mask,     1, mpi_mask, i + 1, 0, MPI_COMM_WORLD );

		MPI_Send( &value,             1, mpi_mask, i + 1, 0, MPI_COMM_WORLD );
				
		//cout << endl << "task # " << mpi_i << " was send to core # " << mpi_i + 1 << endl;
	}

	next_task_index = first_range_tasks_count - 1;
	// send tasks if needed
	// cur_ind_part_control_va and cur_ind_part_work_var are ready from prev for( ; ; ) { } 
	stringstream sstream;

	solved_tasks_count = skip_tasks; // don't solved first skip_tasks problems
	
	// write init info
	WriteSolvingTimeInfo( solving_times, total_solving_times, solved_tasks_count, 
			              sat_count, finding_first_sat_time );
	
	while ( solved_tasks_count < all_tasks_count ) {
		// recieve from core message about solved task 
		// if extra_tasks_count > 0 then this code is unesed
		
		MPI_Recv( &process_sat_count,   1, MPI_INT,          MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status );
		current_status = status;
		MPI_Recv( &solving_times, 1, mpi_solving_time, current_status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status );
		
		solved_tasks_count++;
		next_task_index++;
		if ( verbosity > 0 ) {
			cout << "solved_tasks_count " << solved_tasks_count << endl;
			cout << "next_task_index " << next_task_index << endl;
		}
		
		if ( process_sat_count ) {
			sat_count += process_sat_count;
			cout << "sat_count " << sat_count << endl;
			if ( finding_first_sat_time == 0 ) // first time only
				finding_first_sat_time = MPI_Wtime() - start_time;
		}
		
		WriteSolvingTimeInfo( solving_times, total_solving_times, solved_tasks_count, 
			                  sat_count, finding_first_sat_time );

		if ( process_sat_count && !IsSolveAll )
			break; // exit if SAT set found
		
		if ( next_task_index < all_tasks_count ) {
			// send new index of task
			MPI_Send( &next_task_index, 1, MPI_INT, current_status.MPI_SOURCE, 0, 
				      MPI_COMM_WORLD );
			// send to free core new task in format of minisat input masks

			if ( IsFileAssumptions )
				continue;

			for ( j = 0; j < FULL_MASK_LEN; j++ )
				value[j] = values_arr[next_task_index][j];
				
			MPI_Send( &value, 1, mpi_mask, current_status.MPI_SOURCE, 0, MPI_COMM_WORLD );
		}
	} // while ( solved_tasks_count < all_tasks_count )
	return true;
}

bool MPI_Solver :: ComputeProcessSolve( )
{
	int current_task_index = -1;
	int process_sat_count = 0;
	int current_obj_val = -1;
	MPI_Status status;
	double cnf_time_from_node = 0.0;
	bool IsFirstTaskRecieved = false;
	unsigned int value[FULL_MASK_LEN];
	double solving_times[SOLVING_TIME_LEN];
	minisat22_wrapper m22_wrapper;
	Problem cnf;
	//Solver *S;
	Solver *S;

	for ( unsigned i = 0; i < FULL_MASK_LEN; ++i )
		value[i] = 0;
	extra_tasks_count = 0;

	if ( solver_type == 4 ) { // last version of minisat
		ifstream in( input_cnf_name );
		m22_wrapper.parse_DIMACS_to_problem(in, cnf);
		in.close();
		S = new Solver();
		S->addProblem(cnf);
		S->verbosity        = 0;
		S->IsPredict        = IsPredict;
		S->core_len         = core_len;
		S->start_activity   = start_activity;
		S->max_solving_time = max_solving_time;
		S->max_nof_restarts = max_nof_restarts;
	}
	
	for (;;) {
		// get index of current task
		MPI_Recv( &current_task_index, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status );
		
		if ( current_task_index == -2 )
			MPI_Finalize();

		// with assumptions file we need only current_task_index for reading values from file 
		if ( !IsFileAssumptions ) {
			// first time  recieve from o-rank service information
			if ( !IsFirstTaskRecieved ) {
				// if some CNF extra "doubled" then full_mask = part_mask and they can by
				// different for different values
				MPI_Recv( &extra_tasks_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status );
				if ( extra_tasks_count )
					cout << "\n Recieved extra_tasks_count " << extra_tasks_count << endl;
				if ( !extra_tasks_count ) { // if standart mode do it once
					MPI_Recv( &full_mask, 1, mpi_mask, 0, 0, MPI_COMM_WORLD, &status );
					MPI_Recv( &part_mask, 1, mpi_mask, 0, 0, MPI_COMM_WORLD, &status );
				}
				IsFirstTaskRecieved = true;
			}
			// if extra then take full_mask every time and part_mask == full_mask
			if ( extra_tasks_count ) {
				MPI_Recv( &full_mask, 1, mpi_mask, 0, 0, MPI_COMM_WORLD, &status );
				equalize_arr( part_mask, full_mask );
			}
			// do it anyway
			MPI_Recv( &value, 1, mpi_mask, 0, 0, MPI_COMM_WORLD, &status );
		}

		// Run Solvers
		if ( !SolverRun( S, full_mask, part_mask, value, process_sat_count, current_obj_val,
						 cnf_time_from_node, solving_times, current_task_index ) )
		{ cout << endl << "Error in SolverRun"; return false; }
		
		if ( verbosity > 0 )
			cout << "\n process_sat_count is " << process_sat_count << endl;
		
		MPI_Send( &process_sat_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD );
		MPI_Send( &solving_times, 1, mpi_solving_time, 0, 0, MPI_COMM_WORLD );
	}

	if ( solver_type == 4 )
		delete S;

	return true;
}

//---------------------------------------------------------
bool MPI_Solver :: MPI_ConseqSolve( int argc, char **argv )
{
// read all cnf in current dir and start then in conseq mode
	string dir = string(".");
	int cnf_count = 0;
	vector<string> files;
	vector<string> cnf_files;
	fstream out_file;
	int current_obj_val = -1;
	int process_sat_count= 0;
	double cnf_time_from_node;
	//PBSolver_cut pbs_cut;
	stringstream solve_sstream;
	double start_sec;
	double final_sec;
	Solver *S;

	// MPI start
	//MPI_Request request;
	int corecount = 10, rank = 0;
	/*MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &corecount );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );*/

	files     = vector<string>( );
	cnf_files = vector<string>( );
	// get all files in current dir
	getdir( dir, files );

	for ( unsigned i = 0; i < files.size( ); i++ ) {
		if ( files[i].find( ".cnf" ) != string::npos ) {
			cnf_count++;
			cnf_files.push_back( files[i] );
			if ( rank == 0 ) 
				cout << endl << "founded cnf " << files[i].c_str( );
		}
	}

	if ( cnf_count > corecount ) {
		if ( rank == 0 ) {
			cout << endl << "Warning. Count of cnf-file > corecount";
			cout << endl << "Only first " << corecount << " cnf will be processed";
			cout << endl << "cnf_count changed to corecount";
		}
		cnf_count = corecount;
	}
	else if ( cnf_count == 0 ) {
		if ( rank == 0 ) cout << endl << "Error. No cnf-files in dir";
		return false;
	}

	if ( rank > cnf_count - 1 ) 
		cout << endl << "core # " << rank << " with no job";
	else {// do job
		stringstream sstream;
		sstream << "answer_" << rank + 1;
		string out_file_name = sstream.str( );
		sstream.str( "" );
		sstream.clear();
		
		start_sec = MPI_Wtime( ); // get init time
		input_cnf_name = &cnf_files[rank][0]; // set current name of file
		cout << endl << "input_cnf_name is " << input_cnf_name;
		unsigned int zero_mask[FULL_MASK_LEN];
		for ( int i = 0; i < FULL_MASK_LEN; i++ )
			zero_mask[i] = 0;
		IsHardProblem = 1;
		double solving_times[SOLVING_TIME_LEN];
		if ( !ReadIntCNF( ) ) // Read original CNF
		{ cout << "\n Error in ReadIntCNF" << endl; return 1; }
		cout << endl << "end of ReadIntCNF";
		if ( rank == 0 ) 
			PrintParams( );
		if ( !IsPB ) {
			int current_task_index = 0;
			cout << endl << endl << "Standart mode of SAT solving";
			if ( !SolverRun( S, zero_mask, zero_mask, zero_mask, process_sat_count, 
				             current_obj_val, cnf_time_from_node, solving_times, current_task_index ) )
			{ cout << endl << "Error in SolverRun"; return false; }
			if ( process_sat_count ) {
				if ( !AnalyzeSATset( ) ) {
					// is't needed to deallocate memory - MPI_Abort will do it	
					cout << "\n Error in Analyzer" << endl;
					MPI_Abort( MPI_COMM_WORLD, 0 );
					return false;
				}
			}
		}
		
		final_sec = MPI_Wtime( ) - start_sec;
		
		sstream << input_cnf_name << " " << final_sec << " sec" << endl;
		out_file.open( out_file_name.c_str( ), ios_base :: out ); // open and clear out file
		out_file << sstream.rdbuf( );
		out_file << solve_sstream.rdbuf( );
		cout << endl << "*** sstream " << sstream.str( );
		cout << endl << "*** solve_sstream " << solve_sstream.str( );
		out_file.close( ); 		
	}

	MPI_Finalize( ); // MPI end

	cout << endl << "End of ControlConseqProcessSolve";
	return true;
}

void MPI_Solver :: PrintParams( )
{
	cout << endl << "solver_type is "           << solver_type;              
	cout << endl << "sort_type is "             << sort_type;                
	cout << endl << "koef_val is "              << koef_val;             		
	cout << endl << "schema_type is "		    << schema_type;           
	cout << endl << "full_mask_var_count is "   << full_mask_var_count;  
	cout << endl << "proc_count is "            << corecount;            
	cout << endl << "core_len is "              << core_len;                 
	cout << endl << "start_activity is "        << start_activity;
	cout << endl << "IsConseq is "              << IsConseq;
	cout << endl << "IsHardProblem is "         << IsHardProblem;	         
	cout << endl << "IsPB is "                  << IsPB;                     
	cout << endl << "best_lower_bound is "      << best_lower_bound;        
	cout << endl << "upper_bound is "	        << upper_bound;			  
	cout << endl << "PB_mode is "               << PB_mode;
	cout << endl << "exch_activ is "            << exch_activ;
	cout << endl << "IsSolveAll exch_activ is " << IsSolveAll;
	cout << endl << "max_solving_time "         << max_solving_time;
	cout << endl << "max_nof_restarts "         << max_nof_restarts;
}

//---------------------------------------------------------
bool MPI_Solver :: MPI_Solve( int argc, char **argv )
{
// Solve with MPI
   int temp_corecount = -1,
	   max_possible_tasks_count = 0,
	   process_sat_count = 0,
	   temp_tasks_count = -1,
	   i = 0,
	   mpi_i = 0,
	   first_range_tasks_count = 0;
	unsigned unsigned_one = 1,
				 j = 0;
	unsigned long long int part_var_power = 0;
	double start_sec = 0.0,
		   final_sec = 0.0,
		   whole_time_sec = 0.0;
	// extended variant of full_mask (one more "1")
	unsigned full_mask_ext[FULL_MASK_LEN];

	rank = -1;

	// MPI start
	//MPI_Request request;
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &corecount );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	
	// new type for sending arrays by MPI
	MPI_Type_contiguous( FULL_MASK_LEN, MPI_UNSIGNED, &mpi_mask );
	MPI_Type_commit( &mpi_mask );
	
	// new type with info about solving time of tasks
	MPI_Type_contiguous( SOLVING_TIME_LEN, MPI_DOUBLE, &mpi_solving_time );
	MPI_Type_commit( &mpi_solving_time );

	IsPredict = false;

	if ( corecount < 2 ) 
	{ printf( "\n Error. corecount < 2" ); return false; }

	if ( !ReadIntCNF( ) ) // Read original CNF
	{ cout << "\n Error in ReadIntCNF" << endl; return 1; }

	if ( !MakeVarChoose( ) ) 
	{ cerr << "\n Error in MakeVarChoose" << endl; return false; }

	// get power of 2 that >= corecount
	// 1 core == control core, hence if corecount == 129, temp_corecount == 128
	temp_corecount = 1;
	while ( temp_corecount < ( corecount - 1 ) )
		temp_corecount <<= 1;
	// get maximum possible count of tasks, that >= koef_val*corecount
	max_possible_tasks_count = temp_corecount * koef_val;
	part_mask_var_count = 0;
	// get count of part mask variables
	temp_tasks_count = max_possible_tasks_count;
	while ( temp_tasks_count > 1 ) {
		part_mask_var_count++;
		temp_tasks_count >>= 1;
	}
	if ( part_mask_var_count > var_choose_order.size() )
		part_mask_var_count = var_choose_order.size();
	cout << "part_mask_var_count " << part_mask_var_count << endl;
	// change batch size to treshold value if needed
	if ( var_choose_order.size() - part_mask_var_count > MAX_BATCH_VAR_COUNT ) {
		part_mask_var_count = var_choose_order.size() - MAX_BATCH_VAR_COUNT;
		cout << "part_mask_var_count changed to " << part_mask_var_count << endl;
	}
	if ( part_mask_var_count > MAX_PART_MASK_VAR_COUNT )
		part_mask_var_count = MAX_PART_MASK_VAR_COUNT;
	// get default count of tasks = power of part_mask_var_count
	part_var_power = ( unsigned_one << part_mask_var_count );

	// if default count of tasks < corecount (count of tasks are not enough)
	// and if user want to use extra_tasks
	if ( ( part_var_power < corecount ) && ( part_var_power > ( corecount - 1 ) / 2 ) )
		extra_tasks_count = corecount - ( int )part_var_power - 1;
	else
		extra_tasks_count = 0;

	// for ex if corecount = 8 and part_mask_var_count = 2 then user don't want to use
	// extra_count, or then he would set part_mask_var_count = 3

	all_tasks_count = ( int )( part_var_power ) + extra_tasks_count;

	//cout << "\n***End of ReadIntCNF" << endl;
	//std :: cout << "\n Hello from process # " << rank << endl;
	if ( rank == 0 ) {
		start_sec = MPI_Wtime( ); // get init time

		cout << "*** MPI_Solve is running ***" << endl;
		cout << "max_possible_tasks_count is " << max_possible_tasks_count << endl;
		PrintParams( );
		
		cout << "part_var_power is "    << part_var_power    << endl;
		cout << "extra_tasks_count is " << extra_tasks_count << endl;
		cout << "all_tasks_count is "   << all_tasks_count   << endl;
		if ( skip_tasks >= all_tasks_count ) {
			cerr << "skip_tasks >= all_tasks_count " << endl;
			cerr << skip_tasks << " >= " << all_tasks_count << endl;
			exit;
		}
		
		// part_var_power - count of part variables
		unsigned int **values_arr;
		values_arr = new unsigned int*[( unsigned int )all_tasks_count];
		for( i = 0; i < all_tasks_count; i++ ) { // max length of array is part_var_power
			values_arr[i] = new unsigned int[FULL_MASK_LEN];
			for( j = 0; j < FULL_MASK_LEN; j++ )
				values_arr[i][j] = 0;
		}
		
		if ( !MakeStandartMasks( part_var_power, values_arr ) ) {
			std :: cout << "\n Error in MakeStandartMasks" << endl;
			// is't needed to deallocate memory - MPI_Abort will do it
			MPI_Abort( MPI_COMM_WORLD, 0 );
			return 1;
		}
		cout << "\n Correct end of MakeStandartMasks" << endl;
		
		if ( extra_tasks_count ) { // if more paralleling is needed (ex. for Gifford) 
			orig_tasks_count      = all_tasks_count  - extra_tasks_count;
			full_mask_tasks_count = orig_tasks_count - extra_tasks_count;
			cout << "\n orig_tasks_count is "      << orig_tasks_count      << endl;
			cout << "\n full_mask_tasks_count is " << full_mask_tasks_count << endl;
			if ( !GetExtraTasks( values_arr, full_mask_ext ) ) {
				std :: cout << "\n Error in MakeStandartMasks" << endl;
				MPI_Abort( MPI_COMM_WORLD, 0 );
				return 1;
			}
			cout << "\n GetExtraTasks done";
			cout << "\n full_mask_ext[0] is " << full_mask_ext[0] << endl;
		}
		
		if ( all_tasks_count >= corecount )
			first_range_tasks_count  = corecount - 1;
		else
			first_range_tasks_count  = all_tasks_count;
			//second_range_tasks_count = corecount - first_range_tasks_count - 1;

		cout << "\n first_range_tasks_count is " << first_range_tasks_count << endl;

		if ( !IsPB ) // common CNF mode
			ControlProcessSolve( first_range_tasks_count, full_mask_ext, values_arr );
		// get time of solving
		// int final_sec = cpuTimeInSec( ) - start_sec;
		final_sec = MPI_Wtime( );
		whole_time_sec = final_sec - start_sec;
		cout << endl << "That took %f seconds" << whole_time_sec << endl;
		// write time of solving
		if ( !WriteTimeToFile( whole_time_sec ) ) {
			cout << "\n Error in WriteTimeToFile" << endl;
			MPI_Abort( MPI_COMM_WORLD, 0 );
		}

		for ( i = 0; i < all_tasks_count; i++ )
			delete[] values_arr[i];
		delete[] values_arr;

		cout << "\n Correct deleting of values_arr" << endl;

		// if SAT set was found then write SAT set to file and call MPI_Abort
		if ( sat_count )
			cout << "\n SAT set was found " << endl;
		else
			cout << "\n SAT set was not found " << endl;
		//next_task_index = -1;
		//MPI_Bcast( &next_task_index, 1, MPI_INT, 0, MPI_COMM_WORLD );
			
		// send messages for finalizing
		int break_message = -2;
		for ( int i = 1; i < corecount; i++ )
			MPI_Send( &break_message, 1, MPI_INT, i, 0, MPI_COMM_WORLD );

		MPI_Type_free( &mpi_mask );
		MPI_Type_free( &mpi_solving_time );
		MPI_Finalize( );
		
		//MPI_Abort( MPI_COMM_WORLD, 0 );
	}
	else { // if rank != 0
		if ( !ComputeProcessSolve() ) {
			cout << "\n Error in ComputeProcessSovle" << endl;
			MPI_Abort( MPI_COMM_WORLD, 0 );
			return 1;
		}
	}

	return 0;
}

//---------------------------------------------------------
bool MPI_Solver :: cpuTimeInHours( double full_seconds, int &real_hours, int &real_minutes, int &real_seconds ) 
{
// Time of work in hours, minutes, seconds
	int full_minutes = -1;
	//
	full_minutes = ( int )full_seconds / 60;
	real_seconds = ( int )full_seconds % 60;
	real_hours   = full_minutes / 60;
	real_minutes = full_minutes % 60;
	//
	return true;
}

//---------------------------------------------------------
bool MPI_Solver :: WriteTimeToFile( double whole_time_sec )
{
// Write time of solving to output file
	string line_buffer,
		   TimeIs_char = "\n Time: ",
		   new_line_char = "\n",
		   hour_char = " h ",
		   minute_char = " m ",
		   second_char = " s ";
	int real_hours = -1,
		real_minutes = -1,
		real_seconds = -1;
	
	if ( !cpuTimeInHours( whole_time_sec, real_hours, real_minutes, real_seconds ) ) 
		{ cout << "Error in cpuTimeInHours" << endl; return false; }

	stringstream sstream;
	sstream << TimeIs_char << real_hours << hour_char 
		    << real_minutes << minute_char << real_seconds << second_char;
	
	ofstream solving_info_file;
	solving_info_file.open( solving_info_file_name.c_str( ), ios :: app );
	if ( !solving_info_file.is_open( ) ) {
		std :: cout << "Error in opening of output file " << solving_info_file_name <<  endl;
		return false;
	}
	solving_info_file << sstream.rdbuf( );
	solving_info_file.close( );

	return true;
}