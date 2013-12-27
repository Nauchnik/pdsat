#include "mpi_solver.h"

#pragma warning( disable : 4996 )

const int    MAX_WORD_LENGTH			   = 64;
const int    MAX_LINE_LENGTH               = 524288;
const int    MAX_LINE_LENGTH_2             = 8192;
const int    MEDIUM_STRING_LEN             = 256;
const double TRANSP_COAST                  = 0.000001;
const int    NUM_KEY_BITS                  = 64;

//=============================================================================
// Constructor/Destructor:

MPI_Solver :: MPI_Solver( ) :
	orig_tasks_count            ( 0 ),
	full_mask_tasks_count       ( 0 ),
	exch_activ			        ( 1 ),
	skip_tasks                  ( 0 ),
	base_solving_info_file_name ( "solving_info" ),
	prev_med_time_sum           ( 0 ),
	solving_iteration_count     ( 0 ),
	interrupted_count           ( 0 ),
	max_solving_time_koef       ( 2 ),
	finding_first_sat_time      ( 0 ),
	total_start_time            ( 0 )
{
	solving_times = new double[SOLVING_TIME_LEN];
	for( unsigned i=0; i < SOLVING_TIME_LEN; ++i )
		solving_times[i] = 0;
	total_solving_times.resize( SOLVING_TIME_LEN );
}

MPI_Solver :: ~MPI_Solver( )
{ 
	delete[] solving_times;
}


//---------------------------------------------------------
bool MPI_Solver :: MPI_Solve( int argc, char **argv )
{
// Solve with MPI
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &corecount );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	IsPredict = false;
	double iteration_start_time, iteration_final_time;
	double whole_final_time;
	double whole_start_time = MPI_Wtime( ); // get init time
	stringstream sstream;
	int break_message = -1;
	int stop_message  = -2;
	
	if ( corecount < 2 ) { 
		printf( "Error. corecount < 2" ); MPI_Abort( MPI_COMM_WORLD, 0 );; 
	}
	
	if ( rank != 0 ) {
		if ( !ComputeProcessSolve() ) {
			cerr << "Error in ComputeProcessSovle" << endl; MPI_Abort( MPI_COMM_WORLD, 0 );
		}
	}
	else { // rank == 0
		cout << "*** MPI_Solve is running ***" << endl;
		total_start_time = MPI_Wtime();
		for (;;) {
			sstream << base_solving_info_file_name << "_" << solving_iteration_count;
			solving_info_file_name = sstream.str();
			sstream.clear(); sstream.str("");
			cout << "solving_info_file_name " << solving_info_file_name << endl; 
			
			iteration_start_time = MPI_Wtime();
			ControlProcessSolve();
			iteration_final_time = MPI_Wtime() - iteration_start_time;
			WriteTimeToFile( iteration_final_time );
			solving_iteration_count++;
			max_solving_time *= max_solving_time_koef; // increase time limit

			CollectAssumptionsFiles();
			
			if ( sat_count && !IsSolveAll )
				break; // exit if SAT set found
			
			// send messages for breaking low loop on compute processes
			if ( interrupted_count )
				for ( int i = 1; i < corecount; i++ )
					MPI_Send( &break_message, 1, MPI_INT, i, 0, MPI_COMM_WORLD );
			else break;
		}
		
		whole_final_time = MPI_Wtime( ) - whole_start_time;
		solving_info_file_name = base_solving_info_file_name + "_total";
		WriteTimeToFile( whole_final_time );
		
		// send messages for finalizing
		for ( int i = 1; i < corecount; i++ )
			MPI_Send( &stop_message, 1, MPI_INT, i, 0, MPI_COMM_WORLD );
		
		MPI_Finalize( );
	}

	return 0;
}

void MPI_Solver :: CollectAssumptionsFiles( )
{
// write info from all known assumptions files to 1 new file
	string dir = string(".");
	string str;
	vector<string> files = vector<string>();
	vector<string> :: iterator it;
	getdir( dir, files ); // get all files in current dir
	stringstream sstream;
	sstream << base_known_assumptions_file_name << "_" << (solving_iteration_count-1) << "_rank";
	string old_known_assumptions_file_mask = sstream.str();
	sstream.clear(); sstream.str("");
	sstream << base_known_assumptions_file_name << "_" << solving_iteration_count;
	known_assumptions_file_name = sstream.str();
	sstream.clear(); sstream.str("");
	cout << "known_assumptions_file_name " << known_assumptions_file_name << endl;

	interrupted_count = 0;
	
	ofstream known_assumptions_file;
	known_assumptions_file.open( known_assumptions_file_name.c_str(), ios_base::out );
	if ( !known_assumptions_file.is_open() ) {
		cerr << "!known_assumptions_file.is_open()" << endl;
		cerr << known_assumptions_file_name << endl;
		MPI_Abort( MPI_COMM_WORLD, 0 );
	}
	
	unsigned old_known_assumptions_file_count = 0;
	for ( it = files.begin(); it != files.end(); ++it ) {
		if ( (*it).find( old_known_assumptions_file_mask ) != string::npos ) {
			old_known_assumptions_file_count++;
			ifstream ifile( (*it).c_str() );
			while ( getline( ifile, str ) ) {
				if ( interrupted_count )
					sstream << endl;
				sstream << str;
				interrupted_count++;
			}
			ifile.close();
			known_assumptions_file << sstream.rdbuf();
			sstream.clear(); sstream.str("");
		}
	}
	cout << "old_known_assumptions_file count " << old_known_assumptions_file_count << endl;
	str = "rm " + old_known_assumptions_file_mask + "*";
	system( str.c_str() ); // remove old files
	cout << "system " << str << endl;
	cout << "After CollectAssumptionsFiles() interrupted_count " << interrupted_count << endl;
	
	known_assumptions_file.close();
}

void MPI_Solver :: AddSolvingTimeToArray( ProblemStates cur_problem_state, double cnf_time_from_node, 
										  double *solving_times )
{
	// solving_times[0]  == min
	// solving_times[1]  == max
	// solving_times[2]  == med
	// solving_times[3]  == sat
	if( cur_problem_state != Interrupted ){
		if ( cnf_time_from_node < solving_times[0] )
			solving_times[0] = cnf_time_from_node;
		if ( cnf_time_from_node > solving_times[1] )
			solving_times[1] = cnf_time_from_node;
	}
	switch( cur_problem_state ){ 
		case Solved :
			if ( cnf_time_from_node < 0.0001 ) solving_times[5]++;
			else if ( cnf_time_from_node < 0.001  ) solving_times[6]++;
			else if ( cnf_time_from_node < 0.01 )   solving_times[7]++;
			else if ( cnf_time_from_node < 0.1  )   solving_times[8]++;
			else if ( cnf_time_from_node < 1    )   solving_times[9]++;
			else if ( cnf_time_from_node < 10   )   solving_times[10]++;
			else if ( cnf_time_from_node < 100  )   solving_times[11]++;
			else if ( cnf_time_from_node < 1000 )   solving_times[12]++;
			else solving_times[13]++;
			break;
		case SolvedOnPreprocessing : 
			solving_times[4]++;
			break;
		case Interrupted : 
			solving_times[14]++;
			break;
	}
}

bool MPI_Solver :: SolverRun( Solver *&S, unsigned &local_interrupted_count, int &process_sat_count, 
							  double &cnf_time_from_node, int current_task_index )
{
	if ( verbosity > 1 )
		cout << "start SolverRun()" << endl;

	process_sat_count = 0;
	vec< vec<Lit> > dummy_vec;
	lbool ret;
	double total_time = 0;
	unsigned current_tasks_solved = 0;
	ProblemStates cur_problem_state;
	stringstream sstream, inter_sstream;
	ofstream ofile;
	unsigned batch_interrupted_count = 0;

	sstream << base_known_assumptions_file_name << "_" << solving_iteration_count << "_rank" << rank;
	string new_assumptions_file_name = sstream.str();
	sstream.clear(); sstream.str("");

	solving_times[0] = 1 << 30; // start min len
	for ( unsigned i = 1; i < SOLVING_TIME_LEN; i++ )
		solving_times[i] = 0;
	
	if ( solver_type == 4 ) {
		if ( assumptions_string_count ) // if assumptions in file 
			MakeAssignsFromFile( current_task_index, dummy_vec );
		else
			MakeAssignsFromMasks( full_mask, part_mask, mask_value, dummy_vec );
		
		if ( ( verbosity > 1 ) && ( rank == 1 ) ) {
			cout << "dummy_vec size" << dummy_vec.size() << endl;
			for ( int i = 0; i < dummy_vec.size(); ++i ) {
				for ( int j=0; j < dummy_vec[i].size(); ++j )
					cout << dummy_vec[i][j].x << " ";
				cout << endl;
			}
		}

		for ( int i = 0; i < dummy_vec.size(); ++i )
			if ( dummy_vec[i].size() == 0 ) {
				cerr << "dummy_vec.size() == 0" << endl;
				return false;
			}
		
		uint64_t prev_starts, prev_conflicts, prev_decisions;
		for ( int i=0; i < dummy_vec.size(); ++i ) {
#ifndef _DEBUG
			cnf_time_from_node = MPI_Wtime( );
#endif
			// save current state to check differences
			prev_starts    = S->starts;
			prev_conflicts = S->conflicts;
			prev_decisions = S->decisions;
			
			S->last_time = Minisat :: cpuTime();
			ret = S->solveLimited( dummy_vec[i] );
			
			//ret = S->solveLimited( dummy_vec[i], true, false ); // for SimpSolver
#ifndef _DEBUG
			cnf_time_from_node = MPI_Wtime( ) - cnf_time_from_node;
#endif
			total_time += cnf_time_from_node;

			if ( ret == l_Undef )
				cur_problem_state = Interrupted; // interrupted cause of restarts or time limit
			else if ( ( S->starts - prev_starts <= 1 ) && ( S->conflicts == prev_conflicts ) && ( S->decisions == prev_decisions ) )
				cur_problem_state = SolvedOnPreprocessing;  // solved by BCP
			else
				cur_problem_state = Solved; // just solved
			
			AddSolvingTimeToArray( cur_problem_state, cnf_time_from_node, solving_times );

			if ( cur_problem_state == Interrupted ) {
				batch_interrupted_count++;
				local_interrupted_count++;
				if ( local_interrupted_count > 1 )
					inter_sstream << endl;
				for ( unsigned j = 0; j < dummy_vec[i].size(); ++j )
					inter_sstream << ( dummy_vec[i][j].x % 2 == 0 ) ? "1" : "0";
			}

			if ( ret == l_True ) {
				process_sat_count++;
				cout << "process_sat_count " << process_sat_count << endl;
				if ( !solving_times[3] ) {
					solving_times[3] = cnf_time_from_node; // time of 1st SAT, write only once
					cout << "SAT time " << solving_times[3] << endl;
				}
				b_SAT_set_array.resize( S->model.size() );
				for ( int i=0; i < S->model.size(); i++ )
					b_SAT_set_array[i] = ( S->model[i] == l_True) ? 1 : 0 ;
				// check res file for SAT set existing
				if ( !AnalyzeSATset( ) ) {
					// is't needed to deallocate memory - MPI_Abort will do it	
					cout << "Error in Analyzer" << endl;
					MPI_Abort( MPI_COMM_WORLD, 0 );
					return false;
				}
				if ( !IsSolveAll )
					break;
			}
			//S->loadState();
			//delete S;
			//S->clearDB(); // we don't need to clear DB, because incremental solving is faster
		}
		if ( batch_interrupted_count ) {
			ofile.open( new_assumptions_file_name.c_str(), ios_base::out | ios_base::app );
			ofile << inter_sstream.rdbuf();
			inter_sstream.clear(); inter_sstream.str("");
			ofile.close();
		}
		
		solving_times[2] = total_time / dummy_vec.size(); // med time for current batch
	}
	else 
		{ cout << "solver_type has unknown format" << endl; return false; }

	return true;
}

void MPI_Solver :: WriteSolvingTimeInfo( double *solving_times, unsigned solved_tasks_count )
{
	int time_sec, time_minutes, time_hours;
	
	if ( solving_times[0] < total_solving_times[0] ) // update min time
		total_solving_times[0] = solving_times[0];
	if ( solving_times[1] > total_solving_times[1] ) // update max time
		total_solving_times[1] = solving_times[1];
	prev_med_time_sum += solving_times[2]; // sum of median time for solved batches
	if ( solved_tasks_count > 0 )
		total_solving_times[2] = prev_med_time_sum  / solved_tasks_count; // update median time
	if ( ( total_solving_times[3] == 0 ) && ( solving_times[3] != 0 ) ) // update sat time
		total_solving_times[3] = solving_times[3];
	
	unsigned long long solved_problems_count = 0;
	for ( unsigned i=4; i < SOLVING_TIME_LEN; i++ ) {
		total_solving_times[i] += solving_times[i];
		solved_problems_count  += (unsigned long long)total_solving_times[i];
	}

	stringstream sstream;
	sstream << "assumptions_string_count " << assumptions_string_count << endl;
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
	sstream << "finding_first_sat_time ";
	cpuTimeInHours( finding_first_sat_time, time_hours, time_minutes, time_sec );
	sstream << time_hours << " h " << time_minutes << " m " << time_sec << " s" << endl;
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

bool MPI_Solver :: ControlProcessSolve( )
{
	cout << "ControlProcessSolve is running" << endl;

	if ( solving_iteration_count == 0 ) {
		if ( !ReadIntCNF() ) { // Read original CNF
			cerr << "Error in ReadIntCNF" << endl; MPI_Abort( MPI_COMM_WORLD, 0 );
		}
		if ( !MakeVarChoose() ) { 
			cerr << "Error in MakeVarChoose" << endl; MPI_Abort( MPI_COMM_WORLD, 0 );
		}
	}
	
	ifstream known_assumptions_file( known_assumptions_file_name.c_str() );
	if ( known_assumptions_file.is_open() ) {
		assumptions_string_count = 0;
		cout << "reading of known_assumptions_file_name " << known_assumptions_file_name << endl;
		string str;
		while ( getline( known_assumptions_file, str ) ) {
			if ( str.size() == 0 ) // skip empty strings
				continue;
			if ( str.size() == var_choose_order.size() )
				assumptions_string_count++;
			else {
				cerr << "str.size() != var_choose_order.size()" << endl; 
				cerr << str.size() << " != " << var_choose_order.size() << endl;
				return false;
			}
		}
		cout << "new assumptions_string_count " << assumptions_string_count << endl;
		known_assumptions_file.close();
	}
	else
		cout << "could not open known_assumptions_file " << known_assumptions_file_name << endl;
	
	unsigned max_possible_tasks_count = (unsigned)(pow( 2, ceil( log(corecount - 1)/log(2) ))) * (1<<koef_val);
	cout << "max_possible_tasks_count " << max_possible_tasks_count << endl;
	part_mask_var_count = log(max_possible_tasks_count)/log(2);
	if ( part_mask_var_count > var_choose_order.size() )
		part_mask_var_count = var_choose_order.size();
	
	// change batch size to treshold value if needed
	if ( var_choose_order.size() - part_mask_var_count > RECOMMEND_BATCH_VAR_COUNT ) {
		part_mask_var_count = var_choose_order.size() - RECOMMEND_BATCH_VAR_COUNT;
		cout << "part_mask_var_count changed to " << part_mask_var_count << endl;
	}
	if ( part_mask_var_count > MAX_PART_MASK_VAR_COUNT )
		part_mask_var_count = MAX_PART_MASK_VAR_COUNT;
	if ( var_choose_order.size() - part_mask_var_count > MAX_BATCH_VAR_COUNT ) {
		cerr << "Error. var_choose_order.size() - part_mask_var_count > MAX_BATCH_VAR_COUNT" << endl;
		cerr << var_choose_order.size() - part_mask_var_count << " < " << MAX_BATCH_VAR_COUNT << endl;
		MPI_Abort( MPI_COMM_WORLD, 0 );
	}
	cout << "part_mask_var_count " << part_mask_var_count << endl;

	// get default count of tasks = power of part_mask_var_count
	unsigned part_var_power = ( 1 << part_mask_var_count );
	cout << "part_var_power "    << part_var_power    << endl;
	// TODO add extended tasks counting
	all_tasks_count = part_var_power;
	cout << "all_tasks_count "   << all_tasks_count   << endl;
	
	if ( (int)all_tasks_count < corecount-1 ) {
		cerr << "Error. all_tasks_count < corecount-1" << endl; 
		cerr << all_tasks_count << " < " << corecount-1 << endl;
		MPI_Abort( MPI_COMM_WORLD, 0 );
	}

	if ( ( assumptions_string_count ) && ( assumptions_string_count < all_tasks_count ) ) {
		cout << "solving_iteration_count < all_tasks_count" << endl;
		all_tasks_count = assumptions_string_count;
		cout << "all_tasks_count changed to " << assumptions_string_count << endl;
	}
	
	if ( solving_iteration_count == 0 )
		PrintParams( );
	
	if ( skip_tasks >= all_tasks_count ) {
		cerr << "Error. skip_tasks >= all_tasks_count " << endl;
		cerr << skip_tasks << " >= " << all_tasks_count << endl;
		MPI_Abort( MPI_COMM_WORLD, 0 );
	}
	
	values_arr.resize( all_tasks_count );
	for ( unsigned i = 0; i < values_arr.size(); ++i )
		values_arr[i].resize( FULL_MASK_LEN );
	
	if ( !assumptions_string_count ) {
		if ( !MakeStandardMasks( part_var_power ) ) {
			cerr << "Error in MakeStandartMasks" << endl; MPI_Abort( MPI_COMM_WORLD, 0 );
		}
		cout << "Correct end of MakeStandartMasks" << endl;
	}
	
	unsigned solved_tasks_count = 0;
	// write init info
	WriteSolvingTimeInfo( solving_times, solved_tasks_count );
	
	// send core_len once to every compute process
	for ( int i=0; i < corecount-1; ++i ) {
		MPI_Send( &core_len,                 1, MPI_INT,      i + 1, 0, MPI_COMM_WORLD );
		MPI_Send( &all_tasks_count,          1, MPI_UNSIGNED, i + 1, 0, MPI_COMM_WORLD );
		MPI_Send( &assumptions_string_count, 1, MPI_INT,      i + 1, 0, MPI_COMM_WORLD );
		MPI_Send( &solving_iteration_count,  1, MPI_UNSIGNED, i + 1, 0, MPI_COMM_WORLD );
		MPI_Send( &max_solving_time,         1, MPI_DOUBLE,   i + 1, 0, MPI_COMM_WORLD );
	}
	
	unsigned next_task_index = 0;
	
	unsigned start_tasks_count = min( corecount - 1, (int)all_tasks_count );
	cout << "start_tasks_count " << start_tasks_count << endl;
	// send to all cores (except # 0) tasks from 1st range
	for ( int i = 0; i < start_tasks_count; ++i ) {
		// send new index of task for reading tasks from file
		MPI_Send( &next_task_index, 1, MPI_INT, next_task_index + 1, 0, MPI_COMM_WORLD );
		
		if ( assumptions_string_count == 0 ) { // don't send values when we have file with assimptions
			copy( values_arr[i].begin(), values_arr[i].end(), mask_value );
			MPI_Send( full_mask,  FULL_MASK_LEN, MPI_UNSIGNED, i + 1, 0, MPI_COMM_WORLD );
			MPI_Send( part_mask,  FULL_MASK_LEN, MPI_UNSIGNED, i + 1, 0, MPI_COMM_WORLD );
			MPI_Send( mask_value, FULL_MASK_LEN, MPI_UNSIGNED, i + 1, 0, MPI_COMM_WORLD );
		}
		next_task_index++;
	}

	total_solving_times[0] = 1 << 30; // start min len
	for ( unsigned i = 1; i < total_solving_times.size(); ++i )
		total_solving_times[i] = 0;
	int process_sat_count = 0;
	MPI_Status status, current_status;
	
	while ( solved_tasks_count < all_tasks_count ) {
		// recieve from core message about solved task 		
		MPI_Recv( &process_sat_count, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status );
		current_status = status;
		MPI_Recv( solving_times, SOLVING_TIME_LEN, MPI_DOUBLE, current_status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status );
		solved_tasks_count++;
		if ( verbosity > 0 )
			cout << "solved_tasks_count " << solved_tasks_count << endl;
		
		if ( process_sat_count ) {
			sat_count += process_sat_count;
			cout << "sat_count " << sat_count << endl;
			if ( finding_first_sat_time == 0 ) // first time only
				finding_first_sat_time = MPI_Wtime() - total_start_time;
		}
		
		WriteSolvingTimeInfo( solving_times, solved_tasks_count );
		
		if ( sat_count && !IsSolveAll )
			break; // exit if SAT set found
		
		if ( next_task_index < all_tasks_count ) {
			// send new index of task
			MPI_Send( &next_task_index, 1, MPI_INT, current_status.MPI_SOURCE, 0, MPI_COMM_WORLD );
			if ( assumptions_string_count == 0 ) {
				// send to free core new task in format of minisat input masks
				copy( values_arr[next_task_index].begin(), values_arr[next_task_index].end(), mask_value );
				MPI_Send( mask_value, FULL_MASK_LEN, MPI_UNSIGNED, current_status.MPI_SOURCE, 0, MPI_COMM_WORLD );
			}
			next_task_index++;
		}
	} // while ( solved_tasks_count < all_tasks_count )
	
	return true;
}

bool MPI_Solver :: ComputeProcessSolve( )
{
	minisat22_wrapper m22_wrapper;
	Problem cnf;
	Solver *S;
	MPI_Status status;
	stringstream sstream;
	int current_task_index;
	int process_sat_count = 0;
	double cnf_time_from_node = 0.0;
	bool IsFirstTaskRecieved;
	unsigned local_interrupted_count;
	int val;
	string solving_info_file_name;
	string str;
	ifstream infile;
	
	if ( solver_type == 4 ) { // last version of minisat
		ifstream in( input_cnf_name );
		m22_wrapper.parse_DIMACS_to_problem(in, cnf);
		in.close();
		S = new Solver();
		S->addProblem(cnf);
		S->verbosity        = 0;
		S->IsPredict        = false;
		S->start_activity   = start_activity;
		S->max_nof_restarts = max_nof_restarts;
	}
	
	for (;;) {
		IsFirstTaskRecieved = false;
		assumptions_string_count = 0;
		MPI_Recv( &core_len,                 1, MPI_INT,      0, 0, MPI_COMM_WORLD, &status );
		MPI_Recv( &all_tasks_count,          1, MPI_INT,      0, 0, MPI_COMM_WORLD, &status );
		MPI_Recv( &assumptions_string_count, 1, MPI_INT,      0, 0, MPI_COMM_WORLD, &status );
		MPI_Recv( &solving_iteration_count,  1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, &status );
		MPI_Recv( &max_solving_time,         1, MPI_DOUBLE,   0, 0, MPI_COMM_WORLD, &status );
		local_interrupted_count = 0;
		if ( rank == 1 ) {
			cout << "Received core_len "                 << core_len                 << endl;
			cout << "Received all_tasks_count "          << all_tasks_count          << endl;
			cout << "Received assumptions_string_count " << assumptions_string_count << endl;
			cout << "Received solving_iteration_count "  << solving_iteration_count  << endl;
			cout << "Received max_solving_time "         << max_solving_time         << endl;
		}
		if ( ( assumptions_string_count <= 0 ) && ( solving_iteration_count > 0 ) ) {
			cerr << "Error. ( ( assumptions_string_count <= 0 ) && ( solving_iteration_count > 0 ) )" << endl;
			return false;
		}
		S->core_len         = core_len;
		S->max_solving_time = max_solving_time;
		sstream << base_known_assumptions_file_name << "_" << solving_iteration_count;
		known_assumptions_file_name = sstream.str();
		sstream.clear(); sstream.str("");
		
		if ( solving_iteration_count == 0 ) {
			// read var_choose_order from file
			solving_info_file_name = base_solving_info_file_name + "_0";
			infile.open( solving_info_file_name.c_str() );
			if ( !infile.is_open() ) {
				cerr << "!infile.is_open()" << endl;
				cerr << "solving_info_file_name" << endl;
			}
			bool ReadNext = false;
			while ( getline( infile, str ) ) {
				if ( rank == 1 )
					cout << str << endl;
				if ( str.find( "var_choose_order" ) != string::npos ) {
					ReadNext = true;
					continue;
				}
				if ( ReadNext ) {
					sstream << str;
					if ( rank == 1 )
						cout << "parsing " << str << endl;
					while ( sstream >> val )
						var_choose_order.push_back( val );
					break;
				}
			}
			sstream.clear(); sstream.str("");
			infile.close();
			cout << "rank " << rank << endl;
			cout << "var_choose_order.size() " << var_choose_order.size() << endl;
			for ( unsigned i=0; i < var_choose_order.size(); ++i )
				cout << var_choose_order[i] << " ";
			cout << endl;
		}

		for (;;) {
			// get index of current task
			MPI_Recv( &current_task_index, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status );
			
			if ( current_task_index == -1 )
				break; // stop and get new init values for solving iteration
			else if ( current_task_index == -2 )
				MPI_Finalize( ); // finalize-message from control process
			
			// with assumptions file we need only current_task_index for reading values from file 
			if ( assumptions_string_count == 0 ) {
				if ( !IsFirstTaskRecieved ) {
					MPI_Recv( full_mask, FULL_MASK_LEN, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, &status );
					MPI_Recv( part_mask, FULL_MASK_LEN, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, &status );
					IsFirstTaskRecieved = true;
				}
				MPI_Recv( mask_value, FULL_MASK_LEN, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, &status );
			}

			if ( verbosity > 2 ){
				cout << "full_mask" << endl;
				for( unsigned i=0; i<FULL_MASK_LEN; ++i )
					cout << full_mask[i] << " ";
				cout << endl;
				cout << "part_mask" << endl;
				for( unsigned i=0; i<FULL_MASK_LEN; ++i )
					cout << part_mask[i] << " ";
				cout << endl;
				cout << "mask_value" << endl;
				for( unsigned i=0; i<FULL_MASK_LEN; ++i )
					cout << mask_value[i] << " ";
				cout << endl;
			}
		
			if ( !SolverRun( S, local_interrupted_count, process_sat_count, cnf_time_from_node, current_task_index ) ) { 
				cout << endl << "Error in SolverRun"; return false; 
			}
		
			if ( verbosity > 0 )
				cout << "process_sat_count is " << process_sat_count << endl;
		
			MPI_Send( &process_sat_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD );
			MPI_Send( solving_times, SOLVING_TIME_LEN, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
		}
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
	unsigned local_interrupted_count = 0;

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
		if ( !ReadIntCNF( ) ) // Read original CNF
		{ cout << "\n Error in ReadIntCNF" << endl; return 1; }
		cout << endl << "end of ReadIntCNF";
		if ( rank == 0 ) 
			PrintParams( );
		if ( !IsPB ) {
			int current_task_index = 0;
			cout << endl << endl << "Standart mode of SAT solving";
			if ( !SolverRun( S, local_interrupted_count, process_sat_count, cnf_time_from_node, current_task_index ) )
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
	cout << endl << "koef_val is "              << koef_val;             		
	cout << endl << "schema_type is "		    << schema_type;           
	cout << endl << "full_mask_var_count is "   << full_mask_var_count;  
	cout << endl << "proc_count is "            << corecount;            
	cout << endl << "core_len is "              << core_len;                 
	cout << endl << "start_activity is "        << start_activity;
	cout << endl << "IsConseq is "              << IsConseq;         
	cout << endl << "IsPB is "                  << IsPB;                     
	cout << endl << "best_lower_bound is "      << best_lower_bound;        
	cout << endl << "upper_bound is "	        << upper_bound;			  
	cout << endl << "PB_mode is "               << PB_mode;
	cout << endl << "exch_activ is "            << exch_activ;
	cout << endl << "IsSolveAll exch_activ is " << IsSolveAll;
	cout << endl << "max_solving_time "         << max_solving_time;
	cout << endl << "max_nof_restarts "         << max_nof_restarts;
	cout << endl << "max_solving_time_koef "    << max_solving_time_koef;
	cout << endl;
}

//---------------------------------------------------------
bool MPI_Solver :: WriteTimeToFile( double time_sec )
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

	cout << endl << "That took %f seconds " << time_sec << endl;
	
	cpuTimeInHours( time_sec, real_hours, real_minutes, real_seconds );
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