#include "mpi_predicter.h"

void GetLiteralsFromMasks( unsigned int part_mask[FULL_MASK_LEN], unsigned int value[FULL_MASK_LEN] );

MPI_Predicter :: MPI_Predicter( ) :
	predict_from           ( 0 ),   
	predict_to             ( 0 ), 
	proc_count             ( 1 ),  
	cnf_in_set_count       ( 100 ),
	decomp_set_count       ( 0 ),
	max_var_deep_predict   ( 1 ),
	best_var_num           ( 0 ),
	best_predict_time      ( 0.0 ),
	real_best_var_num      ( 0 ),
	real_best_predict_time ( 0.0 ),
	block_count         ( 0 ),
	slow_cnf_mask_index ( 0 ),
	deep_predict_cur_var( 0 ),
	deep_predict        ( 6 ),
	IsRestartNeeded     ( false ),
	IsDecDecomp         ( false ),
	IsSimulatedGranted  ( false ),
	deep_predict_file_name ( "deep_predict" ),
	var_activity_file_name ( "var_activity" ),
	cur_temperature ( 0 ),
	min_temperature ( 20 ),
	temperature_multiply_koef ( 0.99 ),
	start_temperature_koef ( 0.1 ),
	deep_diff_decomp_set_count ( 100 ),
	point_admission_koef ( 0.2 ),
	global_deep_point_index ( 0 ),
	cur_vars_changing ( 1 ),
	global_checked_points_count ( 0 ),
	global_stopped_points_count ( 0 ),
	global_skipped_points_count ( 0 ),
	IsFirstPoint( true ),
	global_count_var_changing ( 0 ), // init vectors
	predict_file_name( "predict" ),
	record_count ( 0 ),
	ts_strategy ( 0 ),
	max_sat_problems ( 200000000 ),
	decomp_sets_in_block ( 0 ),
	current_predict_start_time ( 0 ),
	current_predict_time ( 0 ),
	whole_get_deep_tasks_time ( 0 ),
	whole_deep_time ( 0 ),
	whole_get_predict_time ( 0 ),
	whole_add_new_unchecked_area_time ( 0 ),
	IsFirstStage ( true ),
	IsRecordUpdated ( false ),
	predict_every_sec ( 2 ),
	max_L2_hamming_distance ( 2 ),
	start_sample_varianse_limit ( 0.000000001 )
{ }

MPI_Predicter :: ~MPI_Predicter( )
{ }

struct points_class
{
    unsigned ones_count;
    unsigned size;
};

bool MPI_Predicter :: ControlProcessPredict( int ProcessListNumber, stringstream &sstream_control )
{
// Predicting of compute cost
	if ( verbosity > 0 ) {
		cout << "Start ControlProcessPredict()" << endl;
		unsigned count = 0;
		for ( unsigned i=0; i < all_tasks_count; ++i )
			if (cnf_start_time_arr[i] > 0)
				count++;
		cout << "count of cnf_start_time_arr[j] > 0 " << count << " from " << all_tasks_count << endl;
	}
	
	sstream_control.str(""); sstream_control.clear();
	sstream_control << "In ControlProcessPredict()" << endl;
	
	if ( verbosity > 0 )
		cout << "all_tasks_count is " << all_tasks_count << endl;
	
	if ( cnf_in_set_count < corecount-1 ) {
		cerr << "Error. cnf_in_set_count < corecount-1" << endl;
		cerr << "too small sample to send first batch of tasks" << endl;
		return false;
	}

	int max_possible_tasks_count = 0,
	   process_sat_count = 0,
	   temp_tasks_count = -1,
	   current_task_index = -1,
	   stop_message = -1,
	   iprobe_message = 0,
	   all_skip_count = 0;
	unsigned int k = 0;
	double whole_time_sec = 0.0,
		   cnf_time_from_node = 0.0;
	int total_sat_count = 0;
	MPI_Status status,
		       current_status;
	MPI_Request request;
	unsigned next_task_index = 0; 
	unsigned part_mask_index = 0;
	solved_tasks_count = 0;

	// send to all cores (except # 0) first tasks and then go to 2nd phase - for every resulted node send new task
	for ( int i = 0; i < corecount-1; i++ ) {	
		// get start part_mask_arr[0]
		if ( next_task_index % cnf_in_set_count == 0 ) { // if new sample then new part_mask
			copy( part_mask_arr[part_mask_index].begin(), part_mask_arr[part_mask_index].end(), part_mask );
			part_mask_index++;
		}
		
		if ( verbosity > 2 ) {
			for ( unsigned j = 0; j < FULL_MASK_LEN; j++ )
				cout << "Sending part_mask"<< j << " " << part_mask[j] << endl;
		}
		
		cnf_start_time_arr[i] = MPI_Wtime( ); // fix current time
		node_list[i] = i + 1; // fix node where SAT problem will be solved
		// send new index of task
		MPI_Send( &i,     1, MPI_INT,  i + 1, ProcessListNumber, MPI_COMM_WORLD );
		MPI_Send( part_mask, FULL_MASK_LEN, MPI_UNSIGNED, i + 1, ProcessListNumber, MPI_COMM_WORLD );
		if ( verbosity > 2 )
			cout << "task # " << i << " was send to core # " << i + 1 << endl;
		next_task_index++;
	}
	if ( verbosity > 0 )
		cout << "Sending of first tasks done" << endl;
	
	double start_sec = MPI_Wtime( ); // get init time
	int int_cur_time = 0, prev_int_cur_time = 0;
	double get_predict_time;
	unsigned cur_get_predict_count = 0;
	// send tasks if needed
	while ( solved_tasks_count < all_tasks_count ) {
		for( ; ; ) { // get predict every PREDICT_EVERY_SEC seconds	
			int_cur_time = ( int )( MPI_Wtime( ) - start_sec );
			if ( ( int_cur_time ) && ( int_cur_time != prev_int_cur_time ) && 
				( int_cur_time % predict_every_sec ) == 0 )
			{
				get_predict_time = MPI_Wtime();
				prev_int_cur_time = int_cur_time;
				if ( !GetPredict( ) ) { 
					cout << "Error in GetPredict " << endl; return false; 
				}

				if ( IsRestartNeeded ) {
					cout << "Fast exit in ControlProcessPredict cause of IsRestartNeeded" << endl;
					IsRestartNeeded = false;
					sstream_control << "IsRestartNeeded true" <<  endl;
					return true;
				}
				
				if ( ( verbosity > 0 ) && ( cnf_to_stop_arr.size() > 0 ) )
					cout << "cnf_to_stop_count " << cnf_to_stop_arr.size() << endl;
				
				// send list of stop-messages
				for ( unsigned i = 0; i < cnf_to_stop_arr.size(); i++ ) {
					MPI_Isend( &stop_message, 1, MPI_INT, node_list[cnf_to_stop_arr[i]], 0, 
							   MPI_COMM_WORLD, &request ); // stop_message == -1
					if ( verbosity > 0 )
						cout << "stop-message was send to node # " 
							 << node_list[cnf_to_stop_arr[i]] << endl;
				}
		
				get_predict_time = MPI_Wtime() - get_predict_time;
				
				while ( ceil(get_predict_time) > predict_every_sec ) {
					predict_every_sec *= 2; // increase treshold  
					cout << "get_predict_time " << get_predict_time << endl;
					cout << "predict_every_sec timed to 2. new value " << predict_every_sec << endl;
				}
				
				if ( cur_get_predict_count++ == 5 ) { // write to file every 5 GetPredict()
					if ( !WritePredictToFile( all_skip_count, whole_time_sec ) ) {
						cout << "Error in WritePredictToFile" << endl;
						MPI_Abort( MPI_COMM_WORLD, 0 );
					}
					cur_get_predict_count = 0;
				}
			}
			MPI_Iprobe( MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &iprobe_message, &status );
			if ( iprobe_message )
				break; // if any message from computing processes, catch it
		} // for( ; ; )
		
		// recieve from core message about solved task    
		MPI_Recv( &current_task_index, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status ); 
		// if 1st message from core # i then get 2nd message from that core
		current_status = status;
		// then get 1 more mes
		MPI_Recv( &process_sat_count,  1, MPI_INT,    current_status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		MPI_Recv( &cnf_time_from_node, 1, MPI_DOUBLE, current_status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		/*double *var_activity = new double[activity_vec_len];
		MPI_Recv( var_activity, activity_vec_len, MPI_DOUBLE, current_status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

		for( unsigned i=0; i < total_var_activity.size(); ++i ) {
			if( ( total_var_activity[i] += var_activity[i] ) > 1e100 )
				for( unsigned j=0; j < total_var_activity.size(); ++j ) // Rescale:
					total_var_activity[j] *= 1e-100;
		}
		delete[] var_activity;*/
		
		// skip old message
		if ( current_status.MPI_TAG < ProcessListNumber ) {
			if ( verbosity > 0 ) {
				cout << "status.MPI_TAG "   << status.MPI_TAG    << endl;
				cout << "ProcessListNumber " << ProcessListNumber << endl;
				cout << "skipped old message" << endl;
			}
			continue;
		}

		// SAT-problem was solved
		cnf_real_time_arr[current_task_index] = cnf_time_from_node; // real time of solving
		
		if ( process_sat_count ) {
			total_sat_count += process_sat_count;
			cout << "total_sat_count " << total_sat_count << endl;
			cnf_status_arr[current_task_index] = 3; // status of CNF is SAT
		}
		else {
			if ( cnf_status_arr[current_task_index] == 0 ) // if status of CNF is not STOPPED
				cnf_status_arr[current_task_index] = 2; // then status of CNF is UNSAT
		}
		
		solved_tasks_count++;
		if ( verbosity > 1 )
			cout << "solved_tasks_count is " << solved_tasks_count << endl;
		
		if ( next_task_index < all_tasks_count ) {
			// skip sending stopped tasks
			while ( ( cnf_status_arr[next_task_index] == 1 ) && ( next_task_index < all_tasks_count ) ) {
				solved_tasks_count++;
				if ( verbosity > 1 ) {
					cout << "skipped sending of task # " << next_task_index    << endl;
					cout << "solved_tasks_count "        << solved_tasks_count << endl;
				}
				next_task_index++;
				//cnf_real_time_arr[next_task_index] = 0.0; // real time of solving
				all_skip_count++;
			}

			// if last tasks were skipped
			if ( next_task_index >= all_tasks_count ) {
				sstream_control << "*next_task_index >= all_tasks_count" << endl;
				sstream_control << next_task_index << " >= " << all_tasks_count << endl;
				break;
			}
			
			cnf_start_time_arr[next_task_index] = MPI_Wtime( ); // set time
			node_list[next_task_index] = current_status.MPI_SOURCE; // get # of node for sending
			
			if ( next_task_index % cnf_in_set_count == 0 ) { // if new sample then new part_mask
				copy( part_mask_arr[part_mask_index].begin(), part_mask_arr[part_mask_index].end(), part_mask );
				part_mask_index++;
			}
			
			if ( verbosity > 2 ) {
				cout << "Sending next_task_index "<< next_task_index << endl;
				for ( unsigned  j = 0; j < FULL_MASK_LEN; j++ )
					cout << "Sending part_mask"<< j << " " << part_mask[j] << endl;
			}
			
			// send new index of task
			MPI_Send( &next_task_index, 1, MPI_INT, current_status.MPI_SOURCE, ProcessListNumber, MPI_COMM_WORLD );
			// send to free core new task in format of minisat input masks
			MPI_Send( part_mask, FULL_MASK_LEN, MPI_UNSIGNED, current_status.MPI_SOURCE, ProcessListNumber, MPI_COMM_WORLD );
			next_task_index++; // if status != STOPPED next_task_index will increase once
		}
	} // while ( solved_tasks_count < tasks_count )
	
	sstream_control << "after while (solved_tasks_count < all_tasks_count) loop" << endl;
	sstream_control << "solved_tasks_count " << solved_tasks_count << endl;
	sstream_control << "global_deep_point_index " << global_deep_point_index << endl;
	sstream_control << "total_decomp_set_count " << global_deep_point_index << endl;
	
	if ( !GetPredict( ) )
		{ cout << "\n Error in GetPredict " << endl; return false;}

	if ( !WritePredictToFile( all_skip_count, whole_time_sec ) ) {
		cout << "\n Error in WritePredictToFile" << endl;
		MPI_Abort( MPI_COMM_WORLD, 0 );
	}

	return true;
}

bool MPI_Predicter :: ComputeProcessPredict( )
{
	// read file with CNF once
	Problem cnf;
	if ( solver_type == 4 ) {
		minisat22_wrapper m22_wrapper;
		ifstream in( input_cnf_name );
		m22_wrapper.parse_DIMACS_to_problem(in, cnf);
		in.close();
	}

	MPI_Status status;
	// get core_len before getting tasks
	MPI_Recv( &core_len,         1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
	MPI_Recv( &activity_vec_len, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
	cout << "Received core_len "         << core_len         << endl;
	cout << "Received activity_vec_len " << activity_vec_len << endl;
	
	int current_task_index;
	double cnf_time_from_node = 0.0;
	int ProcessListNumber;
	vec< vec<Lit> > dummy_vec;
	Solver *S;
	lbool ret;
	unsigned *part_mask_prev = new unsigned[FULL_MASK_LEN];
	double *var_activity = new double[activity_vec_len];
	part_mask_prev[0] = 0; // init to make it differ from part_mask
	bool IsFirstTaskReceived = false;
	bool IsNewPartMask;
	int process_sat_count;

	for (;;) {	
		if ( verbosity > 2 )
			cout << "Before MPI_Recv()" << endl;
		
		do // get index of current task missing stop-messages
		{ 
			if ( verbosity > 0 )
				cout << "on node " << rank << " stop-message was skipped" << endl;
			MPI_Recv( &current_task_index, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		} while ( current_task_index == -1 ); // skip stop-messages

		if ( verbosity > 0 )
			cout << "current_task_index" << current_task_index << endl;
		ProcessListNumber = status.MPI_TAG;
		// receive data from o-rank core in format of minisat input masks
		MPI_Recv( part_mask, FULL_MASK_LEN, MPI_UNSIGNED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		
		if ( verbosity > 2) {
			for ( unsigned i = 0; i < FULL_MASK_LEN; ++i )
				cout << "Received part_mask" << i << " " << part_mask[i] << endl;
		}
		
		IsNewPartMask = false;
		for ( unsigned i = 0; i < FULL_MASK_LEN; ++i )
			if ( part_mask[i] != part_mask_prev[i] ) {
				IsNewPartMask = true;
				break;
			}
		
		if ( IsNewPartMask ) {
			copy( part_mask, part_mask + FULL_MASK_LEN, part_mask_prev );
			if ( solver_type == 4 ) {
				if ( IsFirstTaskReceived ) // if not first time, delete old data
					delete S;
				S = new Solver();
				S->addProblem(cnf);
				S->verbosity        = verbosity;
				S->IsPredict        = IsPredict;
				S->core_len         = core_len;
				S->start_activity   = start_activity;
				S->max_solving_time = max_solving_time;
				S->rank             = rank;
			}
			IsFirstTaskReceived = true;
		}

		// in predict full and part masks are equal
		copy( part_mask, part_mask + FULL_MASK_LEN, full_mask );
		for ( unsigned i = 1; i < FULL_MASK_LEN; i++ )
			mask_value[i] = uint_rand(); // make rand values as assumptions
		
		process_sat_count = 0;
        if ( solver_type == 4 ) {
			MakeAssignsFromMasks( full_mask, part_mask, mask_value, dummy_vec );
			if ( dummy_vec.size() > 1 ) {
				cerr << "Error. In predict mode dummy_vec.size() > 1" << endl;
				cerr << "dummy_vec.size() " << dummy_vec.size() << endl;
				MPI_Finalize( );
			}
			
			cnf_time_from_node = MPI_Wtime( );
			ret = S->solveLimited( dummy_vec[0] );
			cnf_time_from_node = MPI_Wtime( ) - cnf_time_from_node;
			S->GetActivity( var_activity, activity_vec_len ); // get activity of Solver

			if ( cnf_time_from_node < MIN_SOLVE_TIME ) // TODO. maybe 0 - but why?!
				cnf_time_from_node = MIN_SOLVE_TIME;
			if ( ret == l_True ) {
				process_sat_count++;
				cout << "SAT found" << endl;
				cout << "process_sat_count " << process_sat_count << endl;
				b_SAT_set_array.resize( S->model.size() );
				for ( int i=0; i < S->model.size(); i++ )
					b_SAT_set_array[i] = ( S->model[i] == l_True) ? 1 : 0;
				if ( !AnalyzeSATset( ) ) { 	// check res file for SAT set existing
					cout << "Error in Analyzer procedute" << endl;
					return false;
				}
			}
		    S->clearDB();
        }
		else
        { cout << "solver_type has unknown format"; return false; }
		
		MPI_Send( &current_task_index, 1, MPI_INT,    0, ProcessListNumber, MPI_COMM_WORLD );
		MPI_Send( &process_sat_count,  1, MPI_INT,    0, ProcessListNumber, MPI_COMM_WORLD );
		MPI_Send( &cnf_time_from_node, 1, MPI_DOUBLE, 0, ProcessListNumber, MPI_COMM_WORLD );
		//MPI_Send( &var_activity,       1, mpi_var_activity, 0, ProcessListNumber, MPI_COMM_WORLD );
		//MPI_Send( var_activity, activity_vec_len, MPI_DOUBLE, 0, ProcessListNumber, MPI_COMM_WORLD );
	}

	delete[] var_activity;
	delete[] part_mask_prev;
	delete S;
	MPI_Finalize( );
	return true;
}

// ----------------
void GetLiteralsFromMasks( unsigned *part_mask, unsigned *mask_value )
{
// Get onelitetal clauses to make CNF which is hard for solver
	int k = 0;
	stringstream sstream;
	unsigned int cur_var_ind, mask;

	for ( unsigned int i = 1; i < FULL_MASK_LEN; i++ ) {
		for ( int j = 0; j < UINT_LEN; j++ ){
			mask = ( 1 << j );
			if ( part_mask[i] & mask ) {
				cur_var_ind = ( i - 1 ) * UINT_LEN + j + 1;
				if ( mask_value[i] & mask )
					sstream << cur_var_ind;
				else
					sstream << "-" << cur_var_ind;
				sstream << " 0" << endl;
				k++;
			}
		}
	}
	ofstream out_file;
	out_file.open( "literals_slow_cnf", ios :: out );
	out_file << sstream.rdbuf( );
	out_file.close();
}

void MPI_Predicter :: GetInitPoint( )
{
// read point from previous launches
	if ( verbosity > 0 )
		cout << "GetInitPoint() started" << endl;
	fstream deep_predict_file;
	ifstream known_point_file;
	stringstream temp_sstream, sstream;
	known_point_file.open( known_point_file_name.c_str(), ios_base::in );

	if ( known_point_file.is_open() ) { 
		string str;
		while ( getline( known_point_file, str ) ) {
			temp_sstream.str( "" ); temp_sstream.clear( );
			temp_sstream << str;
			temp_sstream >> str; // get first word
			if ( str == "best_predict_time" )
				temp_sstream >> best_predict_time;
			if ( str == "cur_temperature" )
				temp_sstream >> cur_temperature;

			if ( cur_temperature == 0 )
				cur_temperature = best_predict_time * start_temperature_koef;

			if ( best_predict_time > 0.0 )
				IsFirstPoint = false; // don't compute known point

			best_var_num = var_choose_order.size(); 
			// Report about known point once
			sstream << "*** Known init point" << endl << endl;
			sstream << "IsFirstPoint " << IsFirstPoint << endl;
			sstream << "best_predict_time " << best_predict_time << " s" << endl;
			sstream << "cur_temperature " << cur_temperature << endl;
			sstream << "best_var_num " << best_var_num << endl;
			sstream << "best_var_choose_order" << endl;
			for ( unsigned i = 0; i < var_choose_order.size(); i++ )
				sstream << var_choose_order[i] << " ";
			sstream << endl;

			deep_predict_file.open( deep_predict_file_name.c_str(), ios_base::out | ios_base::app );
			deep_predict_file << sstream.rdbuf();
			deep_predict_file.close();
		}

		known_point_file.close();
		/*if ( deep_predict == 3 ) {
		    if ( best_var_num < predict_from ) {
				cout << "Error. best_var_num < predict_from " << endl;
				cout << best_var_num << " < " << predict_from << endl;
				MPI_Abort( MPI_COMM_WORLD, 0 );
		    }
		    if ( best_var_num < predict_to ) {
				cout << "Warning. best_var_num < predict_to " << endl;
				cout << "predict_to changed to best_var_num" << endl;
				predict_to = best_var_num;
		    }
		    else if ( best_var_num > predict_to ) {
				cout << "Warning. best_var_num < predict_to " << endl;
				cout << "predict_to changed to best_var_num" << endl;
				predict_to = best_var_num;
		    }
		}*/
	}
	
	if ( schema_type != "rand" ) { // if schema_type was set by user
		full_mask_var_count = predict_to;
		MakeVarChoose();
	}
	else { // if no file with known point then get random init point
		vector <unsigned> rand_arr;
		if ( core_len < ( unsigned )predict_to ) {
		    core_len = predict_to;
			cout << "core_len changed to predict_to : " << endl;
			cout << core_len << "changed to " << predict_to << endl;
		}
		MakeUniqueRandArr( rand_arr, predict_to, core_len );
		cout << "random init point" << endl;
		var_choose_order.resize( predict_to );
		for ( unsigned i = 0; i < var_choose_order.size(); i++ )
			var_choose_order[i] = ( int )rand_arr[i] + 1;
		
		rand_arr.clear();
	}
	sort( var_choose_order.begin(), var_choose_order.end() );
	cout << "var_choose_order" << endl;
	for ( unsigned i = 0; i < var_choose_order.size(); i++ )
		cout << var_choose_order[i] << " ";
	cout << endl << endl;
}

bool MPI_Predicter :: DeepPredictFindNewUncheckedArea( stringstream &sstream ) 
{
	string str;
	unsigned max_checked, min_hamming_distance;
	bool IsAdding;
	unsigned rand_power_value, rand_L2_start_search_index;
	list<unchecked_area> L2_matches;
	points_class points_class_cur;
	vector<points_class> points_class_vec;
	vector<int> power_values;
	unsigned L2_index;
	unsigned L2_erased_count = 0;
	unsigned L1_erased_count = 0;

	if ( IsRecordUpdated ) {
		sstream << endl << "---Record Updated---" << endl;
		sstream << "solved tasks " << solved_tasks_count << " / " << all_tasks_count << endl;
		sstream << "solved tasks " << (double)( ( solved_tasks_count * 100 ) / all_tasks_count ) << " %" << endl;
		sstream << "cur_vars_changing " << cur_vars_changing << endl;
		sstream << "global_deep_point_index " << global_deep_point_index << endl;
		if ( IsRestartNeeded ) {
			IsRestartNeeded = false;
			sstream << "Warning. IsRestartNeeded == true after end of iteration. Changed to false" << endl;
		}
		IsRecordUpdated = false;
		// set new unchecked area
		boost::dynamic_bitset<> bs = IntVecToBitset( core_len, var_choose_order );
		bool IsCorrectAddingArea = false;
		// remove points from L2 that too far from record
		// and find L2 for current record
		list<unchecked_area> :: iterator L2_it = L2.begin();
		while ( L2_it != L2.end() ) {
			if ( ( (*L2_it).center.count() > bs.count() ) && 
				 ( (*L2_it).center.count() - bs.count() > MAX_DISTANCE_TO_RECORD )
				 )
			{
				L2_it = L2.erase( L2_it );
				L2_erased_count++;
			}
			else {
				if ( (*L2_it).center == bs ) {
					current_unchecked_area = *L2_it;
					IsCorrectAddingArea = true;
				}
				++L2_it;
			}
		}

		list<checked_area> :: iterator L1_it = L1.begin();
		while ( L1_it != L1.end() ) {
			if ( ( (*L1_it).center.count() > bs.count() ) && 
				 ( (*L1_it).center.count() - bs.count() > MAX_DISTANCE_TO_RECORD )
				 )
			{
				L1_it = L1.erase( L1_it );
				L1_erased_count++;
			}
			else ++L1_it;
		}
		
		cout << "L2_erased_count " << L2_erased_count << endl;
		cout << "L1_erased_count " << L1_erased_count << endl;

		if ( !IsCorrectAddingArea )
			sstream << "***Error. !IsCorrectAddingArea" << endl;
		to_string( current_unchecked_area.center, str );
		sstream << "current_unchecked_area center " << endl << str << endl;
		to_string( current_unchecked_area.checked_points, str );
		sstream << "current_unchecked_area checked_points " << endl << str << endl;;
		// make initial values
		cur_vars_changing = 1; // start again from Hamming distance == 1
	} else { // if there were no better points in checked area
		sstream << endl << "---Record not updated---" << endl;
		checked_area c_a;
		list<unchecked_area>::iterator L2_it;
		if ( cur_vars_changing < max_var_deep_predict ) {
			cur_vars_changing++; // try larger Hamming distance
			return true;
		}
		if ( deep_predict == 6 ) { // tabu search mode
			boost::dynamic_bitset<> bs, xor_bs;
			bs = IntVecToBitset( core_len, real_var_choose_order );
			max_checked = 0;
			// find new unchecked area
			sstream << "finding new unchecked_area" << endl;
			sstream << "ts_strategy " << ts_strategy << endl;
			if ( verbosity > 0 )
				cout << "finding new unchecked_area" << endl;
			// find needed criteria and mathces points in neighborhood
			if ( ts_strategy == 0 ) {
				min_hamming_distance = (unsigned)core_len;
				for ( L2_it = L2.begin(); L2_it != L2.end(); L2_it++ ) {
					xor_bs = (*L2_it).center ^ bs;
					if ( xor_bs.count() < min_hamming_distance )
						min_hamming_distance = xor_bs.count(); 
				    IsAdding = true;
				    points_class_cur.ones_count = (*L2_it).center.count();
				    points_class_cur.size = 0;
				    for ( unsigned i=0; i<points_class_vec.size(); i++ ) {
						if ( points_class_cur.ones_count == points_class_vec[i].ones_count ) {
							IsAdding = false;
						    points_class_vec[i].size++;
						    break;    
						}
					}
					if ( IsAdding )
					points_class_vec.push_back( points_class_cur );
				}
				sstream << "min hamming distance from L2 " << min_hamming_distance << endl;
				if ( min_hamming_distance > max_L2_hamming_distance ) {
					cout << "min_hamming_distance > max_L2_hamming_distance " << endl;
					cout << min_hamming_distance << " > " << max_L2_hamming_distance << endl;
					return false;
				}
				// remember mathces
				for ( L2_it = L2.begin(); L2_it != L2.end(); L2_it++ ) {
					xor_bs = (*L2_it).center ^ bs;
					if ( xor_bs.count() == min_hamming_distance )
						L2_matches.push_back( *L2_it );
				}
			}
			ofstream ofile("L2_list", ios_base :: out );
			ofile << "points_class_# : ones_count :  size" << endl;
			for ( unsigned i=0; i < points_class_vec.size(); i++ )
				ofile << i << " : " << points_class_vec[i].ones_count << " : " << points_class_vec[i].size << endl;
			ofile.close();
			points_class_vec.clear();
			if ( L2_matches.size() == 0 ) {
				cout << "Error. L2_matches.size() == 0" << endl;
				exit(1);
			}
			if ( verbosity > 0 )
				cout << "L2_matches.size() " << L2_matches.size() << endl;
			sstream << "L2_matches.size() " << L2_matches.size() << endl;
			if ( verbosity > 0 )
				cout << "L2_matches.size() " << L2_matches.size() << endl;
			for ( L2_it = L2_matches.begin(); L2_it != L2_matches.end(); L2_it++ ) {
				IsAdding = true;	
				for ( unsigned i=0; i < power_values.size(); i++ )
					if ( (*L2_it).center.count() == power_values[i] ) {
						IsAdding = false;
				    	break;    
				    }
				if ( IsAdding )
				power_values.push_back( (*L2_it).center.count() );
			}
			rand_L2_start_search_index = uint_rand() % L2_matches.size();
			rand_power_value           = uint_rand() % power_values.size();
			L2_index = 0;
			IsAdding = false;
			// choose randomly area from randoml class of area center power
			// start search from random part of L2 list, because random_shuffle() don't work for list
			for ( L2_it = L2_matches.begin(); L2_it != L2_matches.end(); L2_it++ ) {
				if ( ( L2_index >= rand_L2_start_search_index ) &&
				     ( (*L2_it).center.count() == power_values[rand_power_value] )
					 )
				{
					IsAdding = true;
					current_unchecked_area = *L2_it;
					break;
				}
				L2_index++;
			}
			// try to find to another side if we hav not found point
			if ( !IsAdding ) {
				L2_index = 0;
				for ( L2_it = L2_matches.begin(); L2_it != L2_matches.end(); L2_it++ ) {
					if ( ( L2_index < rand_L2_start_search_index ) &&
						 ( (*L2_it).center.count() == power_values[rand_power_value] ) ) 
					{
						current_unchecked_area = *L2_it;
						break;
					}
					L2_index++;
				}
			}
			sstream << "power_values ";
			for ( unsigned i=0; i<power_values.size(); i++)
				sstream << power_values[i] << " ";
			sstream << endl;
			sstream << "rand_power_value "           << rand_power_value           << endl;
			sstream << "rand_L2_start_search_index " << rand_L2_start_search_index << endl;
			sstream << "L2_index "					 << L2_index                   << endl;

			power_values.clear();
			L2_matches.clear();
			if ( verbosity > 0 )
				cout << "bofore BitsetToIntVec()" << endl;

			var_choose_order = BitsetToIntVec( current_unchecked_area.center );
			to_string( current_unchecked_area.center, str );
			sstream << "current_unchecked_area center " << endl << str << endl;
			sstream << "current_unchecked_area center len " << endl;
			sstream << current_unchecked_area.center.count() << endl;
			sstream << "var_choose_order " << endl;
			for ( vector<int> :: iterator vec_it = var_choose_order.begin(); 
				vec_it != var_choose_order.end(); vec_it++ )
					sstream << *vec_it << " ";
			sstream << endl;
			if ( current_unchecked_area.center.count() == 0 )
				sstream << "***Error. current_unchecked_area.center.count() == 0" << endl;
		}
	}
	return true;
}

bool MPI_Predicter :: DeepPredictMain( )
{
	// make renadom init point for 1st iteration
	// for other ones as init use best point of prev iteration
	int ProcessListNumber = 0;
	MPI_Request request;
	int stop_message = -1;
	fstream deep_predict_file;
	stringstream sstream;
	whole_deep_time = MPI_Wtime();
	stringstream sstream_control;

	cout << "DeepPredictMain() started" << endl;

	deep_predict_cur_var = predict_to;
	// for equal conditions for every dimension
	global_count_var_changing.resize( max_var_deep_predict );
	for ( unsigned i = 0; i < global_count_var_changing.size(); i++ )
		global_count_var_changing[i] = 0;
	
	// read from file if one exists. if not, get random vector
	GetInitPoint( );
	if ( best_var_num )
		deep_predict_cur_var = best_var_num; // don't use predict_to if known start point
	bool IsFastExit = false;
	double current_time = 0;
	string str;
	
	while (
		   ( ( deep_predict <= 2 ) && 
		     ( deep_predict_cur_var >= predict_from ) 
			 )
			||
		   ( ( deep_predict == 5 ) && 
		     ( ( cur_temperature == 0 ) ||
		       ( cur_temperature > min_temperature ) ) 
			 )
			||
			( ( deep_predict == 3 ) || ( deep_predict == 4 ) || ( deep_predict >= 6 ) )  
		   )	
	{
		if ( IsFirstPoint )
		    current_predict_start_time = MPI_Wtime( );
		
		current_time = MPI_Wtime( );
		GetDeepPredictTasks( );
		if ( verbosity > 0 )
			cout << " GetDeepPredictTasks() done" << endl;
		whole_get_deep_tasks_time += MPI_Wtime( ) - current_time;
		//cout << "GetDeepPredictTasks() time " << MPI_Wtime( ) - current_time << endl;

		current_time = MPI_Wtime( );
		if ( !PrepareForPredict( ) )
			{ cout << "Error in PrepareForPredict" << endl; return false; }
		if ( verbosity > 0 )
			cout << "PrepareForPredict() done" << endl;

		if ( !ControlProcessPredict( ProcessListNumber++, sstream_control ) ) {
			cout << endl << "Error in ControlProcessPredict" << endl;
			MPI_Abort( MPI_COMM_WORLD, 0 );
		}
		if ( verbosity > 0 )
			cout << "ControlProcessPredict() done" << endl;

		deep_predict_file.open( deep_predict_file_name.c_str(), ios_base::out | ios_base::app );
		sstream.str( "" ); sstream.clear( );

		sstream << sstream_control.str();
		
		if ( ( IsFirstPoint ) && ( global_deep_point_index == total_decomp_set_count ) ) {
			sstream << "First point" << endl;
			// set new unchecked area
			current_unchecked_area = *L2.begin();
			to_string( current_unchecked_area.center, str );
			sstream << "current_unchecked_area center " << endl << str << endl;
			to_string( current_unchecked_area.checked_points, str );
			sstream << "current_unchecked_area checked_points " << endl << str << endl;;
			IsFirstPoint = false;
		}
		
		// if method >= 2 and new best point found then repeat search in same dimension
		if ( !DeepPredictFindNewUncheckedArea( sstream ) )
			break;

		if ( var_choose_order.size() == 0 )
			break;
		
		//if ( deep_predict <= 2 ) // end of current iteration if deep_predict <= 2
		//	DecVarLinerDeepMode( sstream, deep_predict_file, mpi_mask, ProcessListNumber );
		global_deep_point_index = 0; // new start of counter for block sending

		deep_predict_file << sstream.rdbuf();
		deep_predict_file.close();
		sstream.str( "" ); sstream.clear( );
		
		// stop all current tasks if not first stage (in this stage all problems must be solve)
		if ( !IsFirstStage ) {
			if ( verbosity > 1 )
				cout << "Extra stop sending " << endl;
			for ( int i = 1; i < corecount; i++ ) {
				MPI_Isend( &stop_message, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &request );
				//cout << "stop-message was send to node # " << i << endl;
			}
		}
		
		if ( IsFastExit )
			break;
	} // while

	if ( ( deep_predict <= 2 ) && ( best_predict_time == real_best_predict_time ) ) {
		best_var_num = real_best_var_num ;
		var_choose_order = real_var_choose_order;
		sstream << "***best point changed to real point " << endl;
	}

	cout << "Deep predict ended" << endl;
	deep_predict_file.open( deep_predict_file_name.c_str(), ios_base::out | ios_base::app );
	
	sort( var_choose_order.begin(), var_choose_order.end() );
	sstream << endl;
	sstream << "best_predict_time "    << best_predict_time << " s" << endl;
	sstream << "best_var_num "         << best_var_num << endl;
	sstream << "best_var_choose_order" << endl;
	for ( int i = 0; i < best_var_num; i++ )
		sstream << real_var_choose_order[i] << " ";
	sstream << endl;
	sstream << "global_count_var_changing" << endl;
	int points_count = 0;
	for ( int i = 0; i < max_var_deep_predict; i++ ) {
		sstream << i + 1 << ":" << global_count_var_changing[i] << " ";
		points_count += global_count_var_changing[i];
	}
	sstream << endl;
	sstream << "points count " << points_count << endl;
	sstream << "temperature " << cur_temperature << endl;
	sstream << "***whole_deep_time "           << MPI_Wtime() - whole_deep_time << " s" << endl;
	sstream << "***whole_get_deep_tasks_time " << whole_get_deep_tasks_time     << " s" << endl;
	sstream << "***whole_get_predict_time "    << whole_get_predict_time        << " s" << endl;
	deep_predict_file << sstream.rdbuf();
	deep_predict_file.close();
	
	global_count_var_changing.clear();
	
	return true;
}

//---------------------------------------------------------
bool MPI_Predicter :: MPI_Predict( int argc, char** argv )
{
// Predicting of computer cost
	stringstream sstream_control;
	rank = -1;

	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &corecount );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	
	IsPredict = true;
	if ( corecount < 2 )
		{ cout << "Error. corecount < 2" << endl; return false; }

	if ( rank == 0 ) { // control node
		cout << "MPI_Predict is running " << endl;
		
		if ( !ReadIntCNF( ) )
		{ cout << "Error in ReadIntCNF" << endl; return 1; }

		if ( predict_to > MAX_CORE_LEN ) {
			cerr << "Warning. predict_to > MAX_PREDICT_TO. Changed to MAX_PREDICT_TO" << endl;
			cerr << predict_to << " > " << MAX_CORE_LEN << endl;
			predict_to = MAX_CORE_LEN;
		}
		if  ( predict_to > (int)core_len ) {
			cout << "core_len changed to predict_to" << endl;
			cout << core_len << " changed to " << predict_to << endl;
			core_len = predict_to;
		}

		if ( ( !predict_to ) && ( core_len ) ) {
			cout << "predict_to changed to core_len" << endl;
			cout << predict_to << " changed to " << core_len << endl;
			predict_to = core_len;
		}

		activity_vec_len = core_len;
		// array of var activity
		total_var_activity.resize( activity_vec_len );
		//for( auto &x : total_var_activity ) x = 0;
		for( vector<double> :: iterator it = total_var_activity.begin(); it != total_var_activity.end(); ++it )
			*it = 0;
		
		// send core_len once to every compute process
		for( int i=0; i < corecount-1; ++i ) {
			MPI_Send( &core_len,         1, MPI_INT,  i + 1, 0, MPI_COMM_WORLD );
			MPI_Send( &activity_vec_len, 1, MPI_INT,  i + 1, 0, MPI_COMM_WORLD );
		}
		
		cout << "verbosity "     << verbosity           << endl;
		cout << "solver_type "   << solver_type         << endl;
		cout << "schema_type "   << schema_type         << endl;
		cout << "corecount "     << corecount           << endl;
		cout << "deep_predict  " << deep_predict        << endl;
		cout << "predict_from "  << predict_from        << endl;
		cout << "predict_to "    << predict_to          << endl;
		cout << "cnf_in_set_count " << cnf_in_set_count << endl; 
		cout << "proc_count "    << proc_count          << endl;
		cout << "core_len "      << core_len            << endl;
		cout << "start_activity "    << start_activity << endl;
		cout << "max_var_deep_predict " << max_var_deep_predict << endl;
		cout << "start_temperature_koef " << start_temperature_koef << endl;
		cout << "point_admission_koef " << point_admission_koef << endl;
		cout << "ts_strategy " << ts_strategy << endl;
		cout << "max_sat_problems " << max_sat_problems << endl;
		cout << "IsFirstStage " << IsFirstStage << endl;
		cout << "max_L2_hamming_distance " << max_L2_hamming_distance << endl;
		
		DeepPredictMain( );
	}
	else { // if rank != 0
		if ( !ComputeProcessPredict( ) ) {
			cout << endl << "Error in ComputeProcessPredict" << endl;
			MPI_Abort( MPI_COMM_WORLD, 0 );
			return false;
		}
	}

	MPI_Abort( MPI_COMM_WORLD, 0 );

	return true;
}

//---------------------------------------------------------
bool MPI_Predicter :: GetRandomValuesArray( unsigned shortcnf_count, vector< vector<unsigned> > &values_arr )
{ 
// get [shortcnf_count] of unique int values that >= 0 and < cnf_count
	// clear after GetIntValuesForMinisat
	for( unsigned i = 0; i < shortcnf_count; i++ )
		for( unsigned j = 0; j < FULL_MASK_LEN; j++ )
			values_arr[i][j] = 0;

	vector< vector<unsigned> > rand_arr;
	unsigned rnd_uint32_count = FULL_MASK_LEN - 1;
	// numbs in rand_arr are long long int, so max 64 vars in decompose set
	MakeRandArr( rand_arr, shortcnf_count, rnd_uint32_count );

	unsigned value_index;
	unsigned mask;
	for( unsigned uint = 0; uint < shortcnf_count; uint++ ) { // for every CNF from random set
		value_index = 0;
		for ( unsigned i = 1; i < FULL_MASK_LEN; i++ ) { // for every 32-bit group of vars
			for( unsigned j = 0; j < UINT_LEN; j++ ) {
				mask = ( 1 << j );
				if ( part_mask[i] & mask ) { // if part_mask bit is 1
					if ( rand_arr[uint][i-1] & mask )
						values_arr[uint][i] += mask;
					value_index++;
				}
			} 
		} 
	} 

	rand_arr.clear();

	for( unsigned uint = 0; uint < shortcnf_count; uint++ ) {
		for( unsigned i = 1; i < FULL_MASK_LEN; i++ ) {
			if ( values_arr[uint][i] )
				values_arr[uint][0] += 1 << ( i-1 );
		}
	}
	
	return true;
}

//---------------------------------------------------------
bool MPI_Predicter :: PrepareForPredict( )
{
// Fill array of tasks for predicting
	unsigned val = 0;

	if ( verbosity > 1 )
		cout << "PrepareForPredict( ) start" << endl;
	
	// array of nodes that solving task with index 0..all_tasks_count-1
	node_list.resize( all_tasks_count );
	// array of CNF status
	cnf_status_arr.resize( all_tasks_count );
	// array of init time of CNF solving
	cnf_start_time_arr.resize( all_tasks_count );
	// array of real time of CNF solving
	cnf_real_time_arr.resize( all_tasks_count );

	for( unsigned i = 0; i < all_tasks_count; i++ ) {
		cnf_status_arr[i]     = 0; // status WAIT
		cnf_start_time_arr[i] = 0.0;
		cnf_real_time_arr[i]  = 0.0;
	}

	// sum times
	sum_time_arr.resize( decomp_set_count );
	// array of median times of CNF in set
	med_time_arr.resize( decomp_set_count );
	// array of predict times
	predict_time_arr.resize( decomp_set_count );
	// array of partly predict times
	predict_part_time_arr.resize( decomp_set_count );
	// array of indexes of array of lengths of sets for modelling
	set_index_arr.resize( decomp_set_count + 1 );
	// array of CNF set status
	set_status_arr.resize( decomp_set_count );
	stopped_cnf_count_arr.resize( decomp_set_count );
	skipped_cnf_count_arr.resize( decomp_set_count );
	solved_cnf_count_arr.resize( decomp_set_count );

	val = 0;
	set_index_arr[0] = 0;
	// set_index_arr = array of CNF set indexes, it depends on set_len_arr
	for( unsigned i = 0; i < decomp_set_count; i++ ) {
		val += set_len_arr[i];
		set_index_arr[i + 1] = val;
		//cout << "\n*** # " << val << endl;
	}

	for( unsigned i = 0; i < decomp_set_count; i++ ) {
		set_status_arr[i]        = 0;
		stopped_cnf_count_arr[i] = 0;
		skipped_cnf_count_arr[i] = 0;
		solved_cnf_count_arr[i]  = 0;
		sum_time_arr[i]          = 0.0;
		med_time_arr[i]          = 0.0;
		predict_time_arr[i]      = 0.0;
		predict_part_time_arr[i] = 0.0;
	}

	return true;
}

void MPI_Predicter :: NewRecordPoint( int set_index )
{
// New record point in deep mode
	fstream deep_predict_file;
	stringstream sstream;
	deep_predict_file.open( deep_predict_file_name.c_str(), ios_base::out | ios_base::app );
	global_count_var_changing[decomp_set_arr[set_index].cur_var_changing - 1]++;
	int old_best_var_num = best_var_num;
	best_var_num = decomp_set_arr[set_index].set_var_count;
	deep_predict_cur_var = best_var_num;
	var_choose_order = decomp_set_arr[set_index].var_choose_order;
	real_var_choose_order = var_choose_order;
	sort( var_choose_order.begin(), var_choose_order.end() );

	double last_predict_record_time = MPI_Wtime() - current_predict_start_time;
	current_predict_time += last_predict_record_time;
	current_predict_start_time = MPI_Wtime( ); // update time

	sstream << "best_predict_time " << best_predict_time << " s" << endl;

	if ( deep_predict == 5 ) { // Simulated annealing
		if ( IsSimulatedGranted ) {
			sstream << "*** Simulated granted" << endl;
			//sstream << "rand_num " << rand_num << endl;
			sstream << "delta " << delta << endl;
			sstream << "exp_value " << exp_value << endl;
			IsSimulatedGranted = false;
		}
		// init or update temperature if == 0 or > limit
		if ( ( cur_temperature == 0 ) || 
		     ( cur_temperature > best_predict_time * start_temperature_koef ) )  
		{     
		    // start temperature
		    cur_temperature = best_predict_time * start_temperature_koef;
		}
		else
			cur_temperature *= temperature_multiply_koef;
		// wtite temperature every new point
		sstream << "cur_temperature " << cur_temperature << endl;
	}

	sstream << "best_var_num " << best_var_num << endl;
	sstream << "best_var_choose_order " << endl;
	for ( int j = 0; j < best_var_num; j++ )
		sstream << var_choose_order[j] << " ";
	sstream << endl;
	sstream << "global_count_var_changing" << endl;
	for ( int i = 0; i < max_var_deep_predict; i++ )
		sstream << i + 1 << ":" << global_count_var_changing[i] << " ";
	sstream << endl << endl;

	if ( IsDecDecomp ) {
		sstream << "real best predict time " << real_best_predict_time << " s" << endl;
		sstream << "real best var num " << real_best_var_num << endl;
		sstream << "real var choose order " << endl;
		for ( int j = 0; j < real_best_var_num; j++ )
			sstream << real_var_choose_order[j] << " ";
		sstream << endl;
	}
	
	if ( 
		( ( deep_predict == 3 ) || ( deep_predict == 5 ) || ( deep_predict == 6 ) )
		&& ( !IsDecDecomp ) && ( !IsFirstPoint ) && ( !IsFirstStage )
	   )
		IsRestartNeeded = true;

	deep_predict_file << sstream.rdbuf();
	deep_predict_file.close();

	sstream.clear(); sstream.str("");
	record_count++;
	sstream << "predict_" << record_count;
	predict_file_name = sstream.str();
	
	if ( !IsFirstStage ) // don't write in first stage - it's expansive
		WritePredictToFile( 0, 0 );
	predict_file_name = "predict";

	ofstream graph_file, var_activity_file;
	if ( IsFirstPoint ) {
		graph_file.open( "graph_file", ios_base :: out ); // erase info from previous launches 
		graph_file << "# best_var_num best_predict_time cnf_in_set_count last_predict_record_time current_predict_time";
		if ( deep_predict == 5 ) // simulated anealing
			graph_file << " cur_temperature";
		graph_file << endl;
		var_activity_file.open( var_activity_file_name.c_str(), ios_base :: out );
	}
	else {
		graph_file.open( "graph_file", ios_base :: app );
		var_activity_file.open( var_activity_file_name.c_str(), ios_base :: app );
	}

	var_activity_file << record_count << endl;
	//for( auto &x : total_var_activity )
	//	var_activity_file << x << " ";
	if ( verbosity > 2 ) {
		cout << "total_var_activity" << endl;
		for( vector<double> :: iterator it = total_var_activity.begin(); it != total_var_activity.end(); ++it ) {
			cout << *it << " ";
			var_activity_file << *it << " ";
		}
		cout << endl;
	}
	
	var_activity_file << endl;
	
	graph_file << record_count << " " << best_var_num << " " << best_predict_time << " " 
		       << cnf_in_set_count << " " << last_predict_record_time << " " << current_predict_time;
	if ( deep_predict == 5 ) 
		graph_file << " " << cur_temperature;

	if ( ( IsFirstStage ) && ( !IsFirstPoint ) && ( best_var_num > old_best_var_num ) ) {
		IsFirstStage = false;
		//cnf_in_set_count *= 10; // increase power of Monte Carlo samples 
		//best_predict_time *= 2; // make new theshold
		cout << "IsFirstStage "      << IsFirstStage      << endl;
		cout << "best_var_num "      << best_var_num      << endl;
		cout << "best_predict_time " << best_predict_time << endl;
		sstream << endl << " *** First stage done ***" << best_predict_time << endl;
		graph_file << " first stage done";
	}

	graph_file << endl;
	graph_file.close();
	var_activity_file.close();
}

bool MPI_Predicter :: IfSimulatedGranted( double predict_time )
{
	if ( deep_predict != 4 )
		return false;
	double rand_num = 0;
	while ( rand_num == 0 ) {
		rand_num = ( unsigned )( uint_rand() % 10 );
		rand_num *= 0.1;
	}
	delta = predict_time - best_predict_time;
	exp_value = exp(-delta / cur_temperature);
	if ( rand_num < exp_value ) {
		IsSimulatedGranted = true;	
		return true;
	}
	else return false;
}

//---------------------------------------------------------
bool MPI_Predicter :: GetPredict()
{
// set_status_arr == 0, if there are some unsolved sat-problem in set 
// set_status_arr == 1, if at least one sat-problem in set was STOPPED
// set_status_arr == 2, if all CNF in set are UNSAT and record
// set_status_arr == 3, if at least one CNF in set is SAT
// set_status_arr == 4, if all CNF in set are UNSAT, not record, but control process couldn't catch to stop it
	int solved_tasks_count = 0;
	unsigned cur_var_num, solved_in_sample_count, 
			 cur_cnf_to_skip_count = 0, 
			 cur_cnf_to_stop_count = 0;
	double cur_predict_time = 0.0,
		   cur_sum_part_time = 0.0,
		   cur_med_part_time = 0.0,
		   tmp_time,
		   time1 = 0,
		   cur_cnf_time = 0;
	unsigned long long temp_llint = 0;
	int set_index_bound = 0,
		cur_cnf_in_set_count = 0;
	double current_time = MPI_Wtime();

	cnf_to_stop_arr.clear(); // every time get stop-list again
	
	// fill arrays of summary and median times in set of CNF
	for ( unsigned i = 0; i < decomp_set_count; i++ ) {
		if ( set_status_arr[i] > 0 ) 
			continue; // skip UNSAT, SAT and STOPPED
		cur_cnf_to_stop_count = 0; // every time create array again
		cur_cnf_to_skip_count = 0;
		sum_time_arr[i]       = 0.0; // init value of sum - every time must start from 0.0
		// count of CNF in set == ( set_index_arr[i + 1] - set_index_arr[i] )
		cur_cnf_in_set_count = set_index_arr[i + 1] - set_index_arr[i];
		solved_in_sample_count = 0;

		tmp_time = MPI_Wtime();
		// if sat-problem is being solved (status 0 or 1) get current time
		for ( unsigned j = set_index_arr[i]; j < set_index_arr[i + 1]; j++ ) {
			switch ( cnf_status_arr[j] ) {
				case 1: // if sat-problem was stopped
					if ( set_status_arr[i] != 3 )
						set_status_arr[i] = 1; // mark status STOP to set
					break;
				case 2:
					solved_in_sample_count++;
					break;
				case 3:
					set_status_arr[i] = 3; // mark status SAT to set
					solved_in_sample_count++;
					break;
				default: // do nothing
					break;
			}
		} // for ( j = set_index_arr[i]; j < set_index_arr[i + 1]; j++ )

		if ( solved_in_sample_count == cur_cnf_in_set_count ) // if all CNF in set has UNSAT status
			set_status_arr[i] = 4; // then mark status UNSAT to set
		
		// further will be only sets whisch are being solved right now
		solved_cnf_count_arr[i] = solved_in_sample_count;

		//max_real_time_sample = 0;
		for ( unsigned j = set_index_arr[i]; j < set_index_arr[i + 1]; j++ ) {
			// if real time from node doesn't exist, use roundly time from 0-core
			if ( cnf_real_time_arr[j] > 0 ) {
				cur_cnf_time = cnf_real_time_arr[j];
				//if ( cur_cnf_time > max_real_time_sample ) 
				//	max_real_time_sample = cur_cnf_time; // compute only real time
			}
			else if ( cnf_start_time_arr[j] > 0 )
				cur_cnf_time = MPI_Wtime( ) - cnf_start_time_arr[j];
			else
				cur_cnf_time = 0; // computing of task wasn't started yet
			// we can stop processing of sample with such task if needed
			sum_time_arr[i] += cur_cnf_time;
		}
		
		if ( deep_predict )
			cur_var_num = decomp_set_arr[i].set_var_count;
		else // if !deep_predict
			cur_var_num = predict_to - i;

		// get current predict time
		med_time_arr[i] = sum_time_arr[i] / (double)cur_cnf_in_set_count;
		cur_predict_time = med_time_arr[i] / (double)proc_count;
		cur_predict_time *= pow( 2, (double)cur_var_num );
		predict_time_arr[i] = cur_predict_time;
		
		// if in set all sat-problems solved and unsat then set can has best predict
		if ( ( set_status_arr[i] == 4  ) && 
			 ( ( best_predict_time == 0.0 ) ||
			   ( ( best_predict_time > 0.0 ) && ( ( predict_time_arr[i] < best_predict_time ) ) ) || 
			   ( ( deep_predict == 5 ) && ( IfSimulatedGranted( predict_time_arr[i] ) ) )
			 )
		   )
		{
			// if all sat-problems solved and unsat then it can be best predict
			set_status_arr[i] = 2;
			IsRecordUpdated = true;
			best_predict_time = predict_time_arr[i];
			
			if ( deep_predict ) // Write info about new point in deep mode
				NewRecordPoint( i );
			//if ( IsRestartNeeded ) // don't stop immidiatly, new record can be found in calculated points
			//	return true;
		}
		// stop, predict >= best
		else if ( ( best_predict_time > 0.0 ) && ( predict_time_arr[i] >= best_predict_time  ) ) {
			if ( ( deep_predict == 5 ) && 
			     ( predict_time_arr[i] < best_predict_time * (1 + point_admission_koef) ) ) // new point can be worst for simalation annealing
			{
				//cout << "passed bad time " << predict_time_arr[i] << " , best time " << best_predict_time << endl; 
				continue;	 
			} 
			
			// don't stop if sample with too simple problems: final - start shows unreal large values
			/*if ( ( max_real_time_sample > 0 ) && ( max_real_time_sample < MIN_STOP_TIME ) ) {
				//cout << "max_time_sample < MIN_STOP_TIME" << endl;
				//cout << max_time_sample << " < " << MIN_STOP_TIME << endl;
				continue;
			}*/
			
			for( unsigned j = set_index_arr[i]; j < set_index_arr[i + 1]; ++j ) {
				if ( ( cnf_start_time_arr[j] > 0 ) && // if solve of sat-problem was started
					 ( cnf_status_arr[j] == 0    ) )  // and still running
				{
					cnf_to_stop_arr.push_back( j ); // then mark it for stopping
					cur_cnf_to_stop_count++;
					cnf_status_arr[j] = 1; // set status STOPPED to CNF
					if ( verbosity > 0 )
						cout << "\n marked STOP to CNF # " << j << endl;
				}
				else if ( cnf_start_time_arr[j] == 0.0 ) { // skip sat-problem if solving wasn't started
					cur_cnf_to_skip_count++;
					cnf_status_arr[j] = 1; // set status STOP to CNF
					//cout << "\n cnf_to_skip_count " << cnf_to_skip_count << endl;
				}
			} // for ( j = set_index_arr[i]; j < set_index_arr[i + 1]; j++ )
			// change stopped and skipped values once
			if ( ( cur_cnf_to_stop_count ) || ( cur_cnf_to_skip_count ) ) {
				set_status_arr[i] = 1; // make status STOP to set of CNF
				if ( cur_cnf_to_stop_count > 0 )
					stopped_cnf_count_arr[i] = cur_cnf_to_stop_count;
				if ( cur_cnf_to_skip_count > 0 )
					skipped_cnf_count_arr[i] = cur_cnf_to_skip_count;
			}

			// if cur_cnf_solved_count == cur_cnf_in_set_count then stopping was late
			// all SAT-problems was solved, but mark set STOP
		}
	} // for ( i = 0; i < decomp_set_count; i++ )
	//cout << "time1 " << time1 << endl;

	stringstream sstream;
	//cout << "before cycle of AddNewUncheckedArea()" << endl;
	if ( deep_predict == 6 ) {
		// add to L2 once for every set
		for ( unsigned i = 0; i < set_status_arr.size(); i++ ) {
			if ( ( decomp_set_arr[i].IsAddedToL2 == false ) && ( set_status_arr[i] > 0 ) )  {
				decomp_set_arr[i].IsAddedToL2 = true;
				if ( set_status_arr[i] == 1 )
					global_stopped_points_count++;
				else 
					global_checked_points_count++;
				
				boost::dynamic_bitset<> bs = IntVecToBitset( core_len, decomp_set_arr[i].var_choose_order );
				AddNewUncheckedArea( bs, sstream );
			}
		}
		fstream deep_predict_file( deep_predict_file_name.c_str(), ios_base::out | ios_base::app );
		deep_predict_file << sstream.rdbuf();
		deep_predict_file.close();
	}

	whole_get_predict_time += MPI_Wtime() - current_time;
	
	return true;
}

//---------------------------------------------------------
bool MPI_Predicter :: WritePredictToFile( int all_skip_count, double whole_time_sec )
{
// Write info about predict to file
	ofstream predict_file;	
	predict_file.open( predict_file_name.c_str( ), ios :: out ); // create and open for writing
	if ( !( predict_file.is_open( ) ) )
	{ cout << "Error in opening of predict_file " << predict_file_name << endl; return false; }

	stringstream sstream;
	sstream << "Predict from "             << predict_from << endl
		    << "Predict to "               << predict_to << endl
			<< "Processor count "          << proc_count << endl
			<< "Count of CNF in set "      << cnf_in_set_count << endl
			<< "Count of CNF in all sets " << all_tasks_count << endl
			<< "Count of core-variables "  << core_len << endl
			<< "Solver type "              << solver_type << endl
			<< "Schema type "              << schema_type << endl
			<< "start activity "           << start_activity << endl;

	if ( deep_predict ) {
		sstream << "\ndeep_predict "     << deep_predict;
		sstream << "\ndecomp_set_count " << decomp_set_count;	
	}
	
	predict_file << sstream.rdbuf( );

	double med_cnf_time, min_cnf_time, max_cnf_time;
	double sample_variance; // sample variance for estimation of derivation
	unsigned max_time_cnf_value_mask = 0;
	for ( unsigned i = 0; i < decomp_set_count; i++ ) {
		med_cnf_time = 0; 
		min_cnf_time = 0; 
		max_cnf_time = 0;
		sample_variance = 0;
		sstream.str( "" );
		sstream.clear();
		sstream << "\n ";

		if ( deep_predict )
			sstream << decomp_set_arr[i].set_var_count;
		else // if !deep_predict
			sstream << predict_to - i;

		sstream << " "  << (unsigned)set_status_arr[i];
		sstream << " sum ";
		sstream.width( 12 );
		sstream.precision( 4 );
		sstream << left << fixed << sum_time_arr[i];
		sstream << " pred ";
		sstream.width( 10 );
		sstream.precision( 2 );
		sstream << left << scientific << predict_time_arr[i];
		
		med_cnf_time = 0;
		// prepare start min and max values
		bool IsFirstNonNullFinded = false;
		unsigned count1 = 0, count2 = 0, count3 = 0, count4 = 0, count5 = 0, count6 = 0, count7 = 0;
		for ( unsigned j = set_index_arr[i]; j < set_index_arr[i + 1]; ++j ) {
			if ( cnf_status_arr[j] <= 1 ) // skip unsolved and stopped
				continue;
			if ( !IsFirstNonNullFinded ) {
				min_cnf_time = max_cnf_time = cnf_real_time_arr[j];
				IsFirstNonNullFinded = true;
			}
			med_cnf_time += cnf_real_time_arr[j];
			if ( cnf_real_time_arr[j] < min_cnf_time )
				min_cnf_time = cnf_real_time_arr[j];
			if ( cnf_real_time_arr[j] > max_cnf_time ) {
				max_cnf_time = cnf_real_time_arr[j];
			}
			if      ( cnf_real_time_arr[j] < 0.01 ) count1++;
			else if ( cnf_real_time_arr[j] < 0.1  ) count2++;
			else if ( cnf_real_time_arr[j] < 1    ) count3++;
			else if ( cnf_real_time_arr[j] < 10   ) count4++;
			else if ( cnf_real_time_arr[j] < 100  ) count5++;
			else if ( cnf_real_time_arr[j] < 1000 ) count6++;
			else count7++;
		}
		
		// compute sample_variance
		if ( solved_cnf_count_arr[i] ) {
			med_cnf_time /= solved_cnf_count_arr[i];
			for ( unsigned j = set_index_arr[i]; j < set_index_arr[i + 1]; j++ ) {
				if ( cnf_status_arr[j] <= 1 )
					continue;
				sample_variance += pow(cnf_real_time_arr[j] - med_cnf_time, 2);
			}
			sample_variance /= ( solved_cnf_count_arr[i] - 1 );
		}

		if ( sample_variance > start_sample_varianse_limit ) {
			cout << "new sample_varianse" << endl;
			cout << "sample_variance " << sample_variance << endl;
			cout << "var_choose_order.size() " << var_choose_order.size() << endl;
			start_sample_varianse_limit *= 2;
		}

		sstream << " min ";
		sstream.width( 15 );
		sstream.precision( 10 );
		sstream << left << fixed << min_cnf_time;
		sstream << " max ";
		sstream.width( 15 );
		sstream.precision( 6 );
		sstream << left << fixed << max_cnf_time;
		sstream << " med ";
		sstream.width( 15 );
		sstream.precision( 10 );
		sstream << left << fixed << med_cnf_time;
		sstream << " s^2 ";
		sstream.width( 15 );
		sstream.precision( 10 );
		sstream << left << fixed << sample_variance;

		sstream << " stopped " << stopped_cnf_count_arr[i];
		sstream << " skipped " << skipped_cnf_count_arr[i];
		sstream << " solved "  << solved_cnf_count_arr[i];
		
		sstream << " (0, 0.01): "   << count1;
		sstream << " (0.01, 0.1): " << count2;
		sstream << " (0.1, 1): "    << count3;
		sstream << " (1, 10): "     << count4;
		sstream << " (10, 100): "   << count5;
		sstream << " (100, 1000): " << count6;
		sstream << " (1000, .): "   << count7;

		predict_file << sstream.rdbuf( );
	}

	sstream.str( "" ); // clear stringstream
	sstream.clear();
	sstream << endl << endl;
	sstream << "All skipped count "       << all_skip_count << endl;
	sstream << "Current solved tasks count " << solved_tasks_count << endl;
	sstream << "Best var num "               << best_var_num;
	sstream << "Best predict time " << left << scientific << best_predict_time << endl;
	sstream << "Predict took time "          << whole_time_sec << endl;
	predict_file << sstream.rdbuf( );
	predict_file.close( );
	
	return true;
}

void MPI_Predicter :: ChangeVarChooseOrder( vector<int> var_choose_order, unsigned change_var_count, 
	                                        unsigned current_var_count, vector<int> &new_var_choose_order )
{
// from given order of vars for choosing make new one with [change_var_count] changed vars
	int non_choose_vars_count = core_len - current_var_count;
	vector<int> non_choose_vars;
	non_choose_vars.resize( non_choose_vars_count );
	
	int k = 0;
	bool IsNewValue;
	// find vars which are not in choose set
	for ( unsigned i = 1; i <= core_len; i++ ) {
		IsNewValue = true;
		for ( unsigned j = 0; j < current_var_count; j++ ) {
			if ( i == var_choose_order[j] ) {
				IsNewValue = false;
				break;
			}
		}
		if ( IsNewValue )
			non_choose_vars[k++] = i;
	}

	vector<unsigned> rand_arr1;
	MakeUniqueRandArr( rand_arr1, change_var_count, current_var_count );
	vector<unsigned> rand_arr2;
	MakeUniqueRandArr( rand_arr2, change_var_count, non_choose_vars_count );

	new_var_choose_order = var_choose_order;
	// add new vars
	for ( unsigned i = 0; i < change_var_count; i++ )
		new_var_choose_order[(unsigned)rand_arr1[i]] = non_choose_vars[(unsigned)rand_arr2[i]];

	sort( new_var_choose_order.begin(), new_var_choose_order.end() );
	
	rand_arr1.clear();
	rand_arr2.clear();
	non_choose_vars.clear();
}

void MPI_Predicter :: GetNewHammingPoint( vector<int> var_choose_order, int cur_vars_changing, int &current_var_count, 
										  vector<int> diff_vec, vector<int> &new_var_choose_order )
{
// make new point with (Hamming distance == change_var_count) from current one 
// differences_vec contains cur_vars_changing '1' in positions where values must be changed
	int cur_var_index = 0;
	new_var_choose_order.clear();
	// if value doesn't exist in diff_vec then add it
	for ( int i = 0; i < current_var_count; i++) {
		if ( find( diff_vec.begin(), diff_vec.end(), var_choose_order[i] ) == diff_vec.end() )
			new_var_choose_order.push_back( var_choose_order[i] );
	}
	for ( unsigned i = 0; i < diff_vec.size(); i++ ) {
		if ( find( var_choose_order.begin(), var_choose_order.end(), diff_vec[i] ) == var_choose_order.end() )
		//if ( !IfValueInArray( diff_vec[i], var_choose_order, current_var_count ) )
			new_var_choose_order.push_back( diff_vec[i] );
	}
	current_var_count = new_var_choose_order.size();
	sort( new_var_choose_order.begin(), new_var_choose_order.end() );
}

void MPI_Predicter :: AllocatePredictArrays( int &cur_tasks_count )
{
	unsigned uint;

	if ( deep_predict_cur_var > MAX_STEP_CNF_IN_SET ) // if too many vars
		cur_tasks_count = cnf_in_set_count;
	else {
		uint = ( 1 << deep_predict_cur_var ); // current count of all tasks
		cur_tasks_count = ( cnf_in_set_count < (int)uint ) ? cnf_in_set_count : uint;
	}
	
	// array of random set lengths 
	set_len_arr.clear();
	for ( unsigned i = 0; i < decomp_set_count; i++ )
	    set_len_arr.push_back( cur_tasks_count );
	
	all_tasks_count = decomp_set_count * cur_tasks_count;
	//part_mask_arr.resize( decomp_set_count );
	part_mask_arr.resize( decomp_set_count );
	//all_values_arr.resize( all_tasks_count );
	
	for( unsigned i = 0; i < part_mask_arr.size(); ++i ) {
		part_mask_arr[i].resize( FULL_MASK_LEN );
		for( unsigned j = 0; j < part_mask_arr[i].size(); ++j )
			part_mask_arr[i][j] = 0;
	}
}

bool MPI_Predicter :: IsPointInCheckedArea( boost::dynamic_bitset<> &point )
{
// check if vector exists in radius of any checked area 
	list< checked_area > :: iterator L1_it;
	boost::dynamic_bitset<> xor_bs;
	for ( L1_it = L1.begin(); L1_it != L1.end(); L1_it++ ) {
		xor_bs = point ^ (*L1_it).center;
		if ( (int)xor_bs.count() <= (*L1_it).radius )
			return true;
	}
	return false;
}

bool MPI_Predicter :: IsPointInUnCheckedArea( boost::dynamic_bitset<> &point )
{
// check if vector exists in radius of any unchecked area 
	list< unchecked_area > :: iterator L2_it;
	boost::dynamic_bitset<> xor_bs;
	unsigned i;
	for ( L2_it = L2.begin(); L2_it != L2.end(); L2_it++ ) {
		xor_bs = point ^ (*L2_it).center;
		if ( (int)xor_bs.none() ) // if new point is center of exist area
			return true;
		else if ( (int)xor_bs.count() == (*L2_it).radius ) {
			//cout << "xor_bs.count() == (*L2_it).radius" << endl;
			// if exists then check with checked points in area
			for ( i = 0; i < xor_bs.size(); i++ )
				if ( xor_bs[i] == 1 ) break;
			// i - index of 1-component
			if ( (*L2_it).checked_points[i] == 1 )
				return true;
		}
	}
	return false;
}

void MPI_Predicter :: AddNewUncheckedArea( boost::dynamic_bitset<> &point, stringstream &sstream )
{
// Add new point as center to list of unchecked areas (L2)
// Add 1 to checked_points vector of new point and all areas from L2 in
// Hamming distanse 1 from new point. Area with current center will be modified too.
// TODO: work with radius > 1
	list<checked_area> :: iterator L1_it;
	list<unchecked_area> :: iterator L2_it;
	boost::dynamic_bitset<> xor_bs;
	unchecked_area new_ua;
	checked_area new_ca;
	unsigned i;
	string str;
	double current_time = MPI_Wtime();

	if ( verbosity > 1 ) {
		cout << "Started AddNewUncheckedArea() with" << endl;
		//sstream << point.to_string() << endl;
	}

	new_ua.center = point;
	new_ua.center_count = new_ua.center.count();
	new_ua.radius = 1;
	new_ua.checked_points.resize( point.size() );

	for ( L1_it = L1.begin(); L1_it != L1.end(); L1_it++ ) {
		// necessary condition - points with close count of 1s
		if ( abs((*L1_it).center_count - new_ua.center_count) <= (*L1_it).radius ) { 
			xor_bs = new_ua.center ^ (*L1_it).center;
			// sufficient condition - points in radius
			if ( (int)xor_bs.count() <= (*L1_it).radius ) {
				for ( i = 0; i < xor_bs.size(); i++ )
					if ( xor_bs[i] == 1 ) break;
				new_ua.checked_points[i] = 1;
			}
		}
	}
	
	if ( verbosity > 1 )
		cout << "In AddNewUncheckedArea() after L1 check" << endl;

	vector< list<unchecked_area> :: iterator > vec_it;
	
	for ( L2_it = L2.begin(); L2_it != L2.end(); L2_it++ ) {
		// necessary condition - points with close count of 1s
		if ( abs((*L2_it).center_count - new_ua.center_count) <= (*L2_it).radius ) { 
			xor_bs = new_ua.center ^ (*L2_it).center;
			// sufficient condition - points in radius
			if ( (int)xor_bs.count() <= (*L2_it).radius ) {
				vec_it.push_back( L2_it );
				for ( i = 0; i < xor_bs.size(); i++ )
					if ( xor_bs[i] == 1 ) break;
				// i - index where 1-component must be set for new area and some areas from L2
				(*L2_it).checked_points[i] = 1;
				new_ua.checked_points[i] = 1;
			}
		}
	}
	
	if ( verbosity > 1 )
		cout << "In AddNewUncheckedArea() after L2 check" << endl;
	
	unsigned move_count = 0;
	// check if some points from L2 must be moved to L1
	// vec_it - iterators to changes points from L2
	for ( unsigned i = 0; i < vec_it.size(); i++ ) {
		L2_it = vec_it[i];
		if ( (*L2_it).checked_points.count() == core_len ) {
			new_ca.center = (*L2_it).center;
			new_ca.center_count = (*L2_it).center.count();
			new_ca.radius = (*L2_it).radius;
			L1.push_back( new_ca );
			if ( current_unchecked_area.center == (*L2_it).center ) {
				sstream << endl << "*** Checked area with Hamming distance " << cur_vars_changing << endl;
				sstream << "best_predict_time " << best_predict_time << endl;
				sstream << endl;
			}
			//to_string( (*L2_it).center, str );
			//sstream << str << endl;
			//sstream << (*L2_it).checked_points.count() << endl;
			L2.erase( L2_it );
			move_count++;
		}
	}
	if ( move_count )
		sstream << "Moved " << move_count << " areas from L2 to L1"  << endl;

	if ( verbosity > 1 )
		cout << "In AddNewUncheckedArea() after L2 -> L1 check" << endl;

	// too expansive to check
	/*bool IsExists = false;
	for ( L2_it = L2.begin(); L2_it != L2.end(); L2_it++ )
		if ( (*L2_it).center == new_ua.center ) {
			sstream << "*** Error. (*L2_it).center == new_ua.center" << endl;
			IsExists = true;
		}

	if ( !IsExists )*/
	L2.push_back( new_ua );
	
	to_string( new_ua.checked_points, str );
	sstream << str << endl;
	if ( ( new_ua.checked_points.count() == 0 ) && ( !IsFirstPoint ) )
		sstream << "***Error. new_ua.center.count() == 0" << endl;
	if ( verbosity > 1 )
		cout << "In AddNewUncheckedArea() end" << endl;
	
	whole_add_new_unchecked_area_time += MPI_Wtime() - current_time ;
}

bool MPI_Predicter :: GetDeepPredictTasks( )
{
// Make tasks for checking neighbors of current point
	vector<int> new_var_choose_order;

	// var_choose_order - current best decomosition (point)
	cout << endl << "*** GetDeepPredictTasks" << endl;

	int points_to_check;
	stringstream sstream;

	cout << "global_deep_point_index " << global_deep_point_index << endl;

	if ( IsFirstPoint )// compute value of function only in first point
		points_to_check = total_decomp_set_count = 1;
	else {
		if ( deep_predict <= 2 )
			decomp_set_count = (cur_vars_changing == 1) ? core_len : max_var_deep_predict * deep_diff_decomp_set_count;
		else { // core_len variants of vectors with Hamming distanse
			// update - make points only for current vars count 
			if ( global_deep_point_index == 0 ) { // calculate only at 1st time
				// get combinations to constuct all vectors on given Hamming distanse
				MakeCombinations( core_len, cur_vars_changing, combinations );
				random_shuffle( combinations.begin(), combinations.end() ); // shuffle vars for getting equal random of changing
				total_decomp_set_count = combinations.size();
				decomp_sets_in_block = max_sat_problems / cnf_in_set_count; // get size of block
				cout << "total_decomp_set_count " << total_decomp_set_count << endl;
				cout << "decomp_sets_in_block " << decomp_sets_in_block << endl;
				if ( decomp_sets_in_block < 1 ) {
					cerr << "decomp_sets_in_block < 1" << endl;
					exit(1);
				}
			}
			// how many points to check
			points_to_check = total_decomp_set_count - global_deep_point_index;
		}
	}

	int cur_tasks_count;
	int cur_var_count; // deep_predict_cur_var may be changed
	unsigned cur_index;
	boost::dynamic_bitset<> point;

	if ( ( deep_predict == 6 ) && ( !IsFirstPoint ) )
		deep_predict_cur_var = current_unchecked_area.center.count();
	
	if ( verbosity > 1 ) {
		cout << "AllocateDeepArrays() done" << endl;
		cout << "decomp_set_count " << decomp_set_count << endl;
	}
	decomp_set_arr.clear();
	decomp_set d_s;
	decomp_set_count = 0;
	int current_skipped = 0;
	//sstream << "current point to check" << endl;
	for ( int i = 0; i < points_to_check; i++ ) { // several decomp sets
		cur_var_count = deep_predict_cur_var; // in Hamming mode (deep_predict >= 3) count can be changed
		if ( IsFirstPoint )
			new_var_choose_order = var_choose_order;
		else {
			cur_index = global_deep_point_index + i;
			if ( deep_predict <= 2 )
				ChangeVarChooseOrder( var_choose_order, cur_vars_changing, cur_var_count, new_var_choose_order );
			else if ( deep_predict != 6 )
				GetNewHammingPoint( var_choose_order, cur_vars_changing, cur_var_count, combinations[cur_index],
									new_var_choose_order );

			if ( deep_predict == 6 ) {
				if ( current_unchecked_area.checked_points[cur_index] == 1 ) {
					current_skipped++;
					continue; // checked already
				}
				else {
					point = current_unchecked_area.center;
					point[cur_index] = ( current_unchecked_area.center[cur_index] == 1 ) ? 0 : 1;
				}
				//sstream << point.to_string() << endl;
				new_var_choose_order = BitsetToIntVec( point );
				cur_var_count = new_var_choose_order.size();
			}
		}
		
		d_s.var_choose_order = new_var_choose_order;
		//PrintVector( d_s.var_choose_order );
		d_s.cur_var_changing = cur_vars_changing;
		d_s.set_var_count = cur_var_count;
		d_s.IsAddedToL2 = false;
		decomp_set_arr.push_back( d_s ); // add new decomp set
		
		decomp_set_count++;
		if ( decomp_set_count == decomp_sets_in_block )
			break;
	} // for ( int i = 0; i < decomp_set_count; i++ )
	random_shuffle( decomp_set_arr.begin(), decomp_set_arr.end() ); // shuffle to erase not noly first variables in decomp set

	if ( verbosity > 0 )
		cout << "after decomp_set_arr creating" << endl;

	global_skipped_points_count += current_skipped;

	// Common for all procedures of getting deep tasks
	AllocatePredictArrays( cur_tasks_count );

	if ( verbosity > 0 )
		cout << "AllocatePredictArrays() done " << endl;
	
	sstream << endl << "***deep_predict_cur_var " << deep_predict_cur_var << endl;
	sstream << "total_decomp_set_count " << total_decomp_set_count << endl;
	sstream << "global_deep_point_index " << global_deep_point_index << endl;
	sstream << "decomp_sets_in_block " << decomp_sets_in_block << endl;
	sstream << "must be checked " << decomp_set_count << " from " << points_to_check << endl;
	//sstream << "points_to_check " << points_to_check << endl;
	//sstream << "decomp_set_count " << decomp_set_count << endl;
	sstream << "current_skipped " << current_skipped << endl;
	//sstream << "cur_vars_changing " << cur_vars_changing << endl;
	sstream << "all_tasks_count " << all_tasks_count << endl;
	sstream << "L1.size() " << L1.size() << endl;
	sstream << "L2.size() " << L2.size() << endl;
	sstream << "global_checked_points_count " << global_checked_points_count << endl;
	sstream << "global_stopped_points_count " << global_stopped_points_count << endl;
	sstream << "global_skipped_points_count " << global_skipped_points_count << endl;
	sstream << "whole_deep_time "             << MPI_Wtime() - whole_deep_time << " s" << endl;
	sstream << "whole_get_deep_tasks_time "   << whole_get_deep_tasks_time     << " s" << endl;
	sstream << "whole_get_predict_time "      << whole_get_predict_time       << " s" << endl;
	sstream << "whole_add_new_unchecked_area_time " << whole_add_new_unchecked_area_time << " s" << endl;
	
	sstream << "best_var " << real_var_choose_order.size() << endl;
	if ( real_var_choose_order != var_choose_order ) {
		sstream << "* real_var_choose_order " << endl;
		for ( unsigned i = 0; i < real_var_choose_order.size(); i++ )
			sstream << real_var_choose_order[i] << " ";
		sstream << endl;
		sstream << "best_predict_time " << best_predict_time << endl;
	}
	else {
		sstream << "var_choose_order " << endl;
		for ( unsigned i = 0; i < var_choose_order.size(); i++ )
			sstream << var_choose_order[i] << " ";
		sstream << endl;
		sstream << "best_var " << var_choose_order.size() << endl;
		sstream << "best_predict_time " << best_predict_time << endl;
	}

	int all_values_arr_index = 0;

	if ( verbosity > 0 )
		cout << sstream.str();

	// Make predict part and value mask arrays
	for ( unsigned i = 0; i < decomp_set_arr.size(); ++i ) {
		// get current [part_mask]
		part_mask_var_count = full_mask_var_count = decomp_set_arr[i].set_var_count;
		if ( !GetMainMasksFromVarChoose( decomp_set_arr[i].var_choose_order ) )
			{ cout << "Error in GetMainMasksFromVarChoose" << endl; return false; }
		for ( unsigned j = 0; j < FULL_MASK_LEN; j++ )
			part_mask_arr[all_values_arr_index][j] = part_mask[j];
		all_values_arr_index++;
	}
	
	global_deep_point_index += decomp_set_count + current_skipped;
	
	if ( verbosity > 0 )
		cout << "After main loop in GetDeepPredictTasks()" << endl;
	
	sstream << endl;
	fstream deep_predict_file( deep_predict_file_name.c_str(), ios_base::out | ios_base::app );
	deep_predict_file << sstream.rdbuf();
	deep_predict_file.close();

	return true;
}