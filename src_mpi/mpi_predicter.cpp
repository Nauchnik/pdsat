#include "mpi_predicter.h"
#include "minisat/utils/SolverStateAccessor.h"

MPI_Predicter :: MPI_Predicter( ) :
	predict_from           ( 0 ),   
	predict_to             ( 0 ), 
	proc_count             ( 1 ),  
	max_var_deep_predict   ( 1 ),
	best_var_num           ( 0 ),
	best_predict_time      ( 0.0 ),
	real_best_var_num      ( 0 ),
	real_best_predict_time ( 0.0 ),
	block_count         ( 0 ),
	deep_predict        ( 6 ),
	IsRestartNeeded     ( false ),
	IsDecDecomp         ( false ),
	isSimulatedGranted  ( false ),
	deep_predict_file_name ( "deep_predict" ),
	var_activity_file_name ( "var_activity" ),
	cur_temperature ( 0 ),
	min_temperature ( 20 ),
	temperature_multiply_koef ( 0.9 ),
	start_temperature_koef ( 0.1 ),
	deep_diff_decomp_set_count ( 100 ),
	point_admission_koef ( 0.2 ),
	global_deep_point_index ( 0 ),
	cur_vars_changing ( 1 ),
	global_checked_points_count ( 0 ),
	global_stopped_points_count ( 0 ),
	global_skipped_points_count ( 0 ),
	isFirstPoint( true ),
	global_count_var_changing ( 0 ), // init vectors
	predict_file_name( "predict" ),
	record_count ( 0 ),
	ts_strategy ( 0 ),
	current_predict_start_time ( 0 ),
	current_predict_time ( 0 ),
	whole_get_deep_tasks_time ( 0 ),
	whole_deep_time ( 0 ),
	whole_get_predict_time ( 0 ),
	whole_add_new_unchecked_area_time ( 0 ),
	IsFirstStage ( false ),
	IsRecordUpdated ( false ),
	predict_every_sec ( 2 ),
	max_L2_hamming_distance ( 2 ),
	start_sample_variance_limit ( 0.000000001 ),
	evaluation_type ( "time" ),
	best_cnf_in_set_count ( 0 ),
	unupdated_count ( 0 ),
	prev_area_best_predict_time ( 0 ),
	blob_var_count ( 0 ),
	cur_point_number ( 0 ),
	tmp_cnf_process_name ( "" ),
	isSolvedOnPreprocessing ( 0 ),
	best_solved_in_time ( 0 ),
	best_time_limit ( 0.0 ),
	predict_time_limit_step ( 0 ),
	template_cnf_size ( 0 )
{
	array_message = NULL;
	for( unsigned i=0; i < PREDICT_TIMES_COUNT; i++ )
		best_predict_time_arr[i] = 0.0;
}

MPI_Predicter :: ~MPI_Predicter( )
{ }

// comparison for sorting
bool ua_compareByActivity(const unchecked_area &a, const unchecked_area &b)
{
	return a.med_var_activity > b.med_var_activity;
}

bool ds_compareByActivity(const decomp_set &a, const decomp_set &b)
{
	return a.med_var_activity > b.med_var_activity;
}

void MPI_Predicter :: SendPredictTask( int ProcessListNumber, int process_number_to_send, int &cur_task_index, unsigned &cur_decomp_set_index )
{
	if ( verbosity > 1 )
		cout << "SendPredictTask() start" << endl;
	
	cnf_start_time_arr[cur_task_index] = Minisat::cpuTime(); // fix current time
	node_list[cur_task_index] = process_number_to_send; // fix node where SAT problem will be solved
	
	if ( ( !cur_task_index ) || ( cur_task_index % (int)cnf_in_set_count == 0 ) ) { // if new sample then new set to all precesses
		if ( verbosity > 1 ) {
			cout << "In SendPredictTask() new decomp set" << endl; 
			cout << "cur_task_index " << cur_task_index << " cnf_in_set_count " << cnf_in_set_count << endl;
		}
		if ( cur_decomp_set_index > decomp_set_arr.size() - 1 ) {
			cerr << "cur_decomp_set_index > decomp_set_arr.size() - 1" << endl;
			cerr << cur_decomp_set_index << " > " << decomp_set_arr.size() - 1 << endl;
			MPI_Abort( MPI_COMM_WORLD, 0 );
		}
		if ( array_message )
			delete[] array_message;
		array_message_size = decomp_set_arr[cur_decomp_set_index].var_choose_order.size() + 1; // first cell for index
		array_message = new int[array_message_size];
		for ( unsigned i=1; i < array_message_size; ++i ) 
			array_message[i] = decomp_set_arr[cur_decomp_set_index].var_choose_order[i-1];
		for ( unsigned i=0; i < vec_IsDecompSetSendedToProcess.size(); ++i ) 
			vec_IsDecompSetSendedToProcess[i] = false;
		if ( verbosity > 1 ) {
			cout << "Ready to send array_message with size " << array_message_size << endl;
		}
		cur_decomp_set_index++;
	}
	
	// send decomp_set if needed
	if ( vec_IsDecompSetSendedToProcess[process_number_to_send-1] == false ) {
		array_message[0] = cur_task_index;
		MPI_Send( array_message, array_message_size, MPI_INT, process_number_to_send, ProcessListNumber, MPI_COMM_WORLD );
		vec_IsDecompSetSendedToProcess[process_number_to_send-1] = true;
		if ( verbosity > 1 )
			cout << "Sending array_message with size " << array_message_size << endl;
	}
	else {
		MPI_Send( &cur_task_index, 1, MPI_INT, process_number_to_send, ProcessListNumber, MPI_COMM_WORLD );
		if ( verbosity > 1 )
			cout << "Sending cur_task_index " << cur_task_index << endl;
	}
	
	if ( verbosity > 1 )
		cout << "SendPredictTask() done" << endl;
	
	cur_task_index++;
}

bool MPI_Predicter :: ControlProcessPredict( int ProcessListNumber, stringstream &sstream_control )
{
// Predicting of compute cost
	if ( verbosity > 0 ) {
		cout << "Start ControlProcessPredict()" << endl;
		unsigned count = 0;
		for ( unsigned i=0; i < all_tasks_count; ++i )
			if ( cnf_start_time_arr[i] > 0 ) count++;
		cout << "count of cnf_start_time_arr[j] > 0 " << count << " from " << all_tasks_count << endl;
	}
	
	sstream_control.str(""); sstream_control.clear();
	sstream_control << "In ControlProcessPredict()" << endl;
	
	if ( te > 0.0 ) {
		predict_time_limites.resize( 100 );
		predict_time_limit_step = te / (double)predict_time_limites.size();
		predict_time_limites[predict_time_limites.size() - 1] = te;
		for ( int i = predict_time_limites.size() - 2; i >= 0; --i )
			predict_time_limites[i] = predict_time_limites[i+1] - predict_time_limit_step;
	}
	
	if ( verbosity > 0 )
		std::cout << "all_tasks_count is " << all_tasks_count << std::endl;
	
	if ( (int)cnf_in_set_count < corecount-1 ) {
		std::cerr << "Error. cnf_in_set_count < corecount-1" << std::endl;
		std::cerr << "too small sample to send first batch of tasks" << std::endl;
		MPI_Abort( MPI_COMM_WORLD, 0 );
	}
	
	int cur_task_index = 0; 
	unsigned cur_decomp_set_index = 0;
	unsigned part_mask_index = 0;
	
	// send to all cores (except # 0) first tasks and then go to 2nd phase - for every resulted node send new task
	for ( int i = 0; i < corecount-1; ++i )
		SendPredictTask( ProcessListNumber, i+1, cur_task_index, cur_decomp_set_index );
	
	std::cout << "Sending of first tasks done" << std::endl;
	
	int task_index_from_node = -1,
	   stop_message = -1,
	   iprobe_message = 0,
	   all_skip_count = 0;
	double whole_time_sec = 0.0;
	int total_sat_count = 0;
	MPI_Status status,
		       current_status;
	MPI_Request request;
	solved_tasks_count = 0;
	cnf_time_from_node = 0.0;
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
				get_predict_time = Minisat::cpuTime();
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
				
				if ( !isSolverSystemCalling ) {
					if ( ( verbosity > 0 ) && ( cnf_to_stop_arr.size() > 0 ) )
						std::cout << "cnf_to_stop_count " << cnf_to_stop_arr.size() << std::endl;
					// send list of stop-messages
					for ( unsigned i = 0; i < cnf_to_stop_arr.size(); i++ ) {
						MPI_Isend( &stop_message, 1, MPI_INT, node_list[cnf_to_stop_arr[i]], 0, 
								   MPI_COMM_WORLD, &request ); // stop_message == -1
						if ( verbosity > 0 )
							std::cout << "stop-message was send to node # " 
								      << node_list[cnf_to_stop_arr[i]] << std::endl;
					}
				}
				
				get_predict_time = Minisat::cpuTime() - get_predict_time;
				
				while ( ceil(get_predict_time) > predict_every_sec ) {
					predict_every_sec *= 2; // increase treshold  
					std::cout << "get_predict_time " << get_predict_time << std::endl;
					std::cout << "predict_every_sec timed to 2. new value " << predict_every_sec << std::endl;
				}
				
				if ( cur_get_predict_count++ == 5 ) { // write to file every 5 GetPredict()
					if ( !WritePredictToFile( all_skip_count, whole_time_sec ) ) {
						std::cout << "Error in WritePredictToFile" << std::endl;
						MPI_Abort( MPI_COMM_WORLD, 0 );
					}
					cur_get_predict_count = 0;
				}
			}
			MPI_Iprobe( MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &iprobe_message, &status );
			if ( iprobe_message )
				break; // if any message from computing processes, catch it
		} // for( ; ; )
		
		if ( verbosity > 1 )
			std::cout << "In ControlProcess() before recieving results" << std::endl;
		
		// recieve from core message about solved task    
		MPI_Recv( &task_index_from_node,       1, MPI_INT,    MPI_ANY_SOURCE,            MPI_ANY_TAG, MPI_COMM_WORLD, &status ); 
		//cout << "recv current_task_index " << current_task_index << endl;
		// if 1st message from core # i then get 2nd message from that core
		current_status = status;
		// then get 1 more mes
		MPI_Recv( &process_sat_count,          1, MPI_INT,    current_status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		MPI_Recv( &cnf_time_from_node,         1, MPI_DOUBLE, current_status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		if ( !isSolverSystemCalling ) {
			MPI_Recv( &isSolvedOnPreprocessing,    1, MPI_INT, current_status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
			MPI_Recv( var_activity, activity_vec_len, MPI_DOUBLE, current_status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
			for( unsigned i=0; i < total_var_activity.size(); ++i ) {
				if( ( total_var_activity[i] += var_activity[i] ) > 1e50 )
					for( unsigned j=0; j < total_var_activity.size(); ++j ) // Rescale:
						total_var_activity[j] *= 1e-50;
			}
		}
		
		if ( verbosity > 1 )
			cout << "Received result with current_task_index " << task_index_from_node << endl;
		
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
		cnf_real_time_arr[task_index_from_node] = cnf_time_from_node; // real time of solving
		cnf_prepr_arr[task_index_from_node] = isSolvedOnPreprocessing;
		
		if ( process_sat_count == 1 ) {
			total_sat_count += process_sat_count;
			//cout << "total_sat_count " << total_sat_count << endl;
			cnf_status_arr[task_index_from_node] = 3; // status of CNF is SAT
			cnf_issat_arr[task_index_from_node] = true;
		}
		else if ( process_sat_count == 0 ) {
			if ( cnf_status_arr[task_index_from_node] == 0 ) // if status of CNF is not STOPPED
				cnf_status_arr[task_index_from_node] = 2; // then status of CNF is UNSAT
		}
		else if ( process_sat_count > 1 ) {
			cerr << "process_sat_count > 1" << endl;
			exit(1);
		}
		
		solved_tasks_count++;
		if ( verbosity > 1 ) {
			cout << "solved_tasks_count is " << solved_tasks_count << endl;
			cout << "cur_task_index " << cur_task_index << endl;
		}

		/*unsigned sat_count = 0;
		for ( auto x = cnf_issat_arr.begin(); x != cnf_issat_arr.end(); x++ )
			if (*x) sat_count++;
		cout << endl << "solved_tasks_count " << solved_tasks_count << endl;
		cout << "sat_count " << sat_count << endl;
		cout << "task_index_from_node " << task_index_from_node << endl;
		cout << "cnf_time_from_node " << cnf_time_from_node << endl; 
		cout << "current_status.MPI_SOURCE " << current_status.MPI_SOURCE << endl;*/
		
		if ( (unsigned)cur_task_index < all_tasks_count ) {
			// skip sending stopped tasks
			while ( ( cnf_status_arr[cur_task_index] == 1 ) && ( (unsigned)cur_task_index < all_tasks_count ) ) {
				solved_tasks_count++;
				all_skip_count++;
				cur_task_index++;
				if ( verbosity > 1 ) {
					cout << "skipped sending of task # " << cur_task_index    << endl;
					cout << "solved_tasks_count "        << solved_tasks_count << endl;
				}
			}
			
			// if last tasks were skipped
			if ( (unsigned)cur_task_index >= all_tasks_count ) {
				sstream_control << "cur_task_index >= all_tasks_count" << endl;
				sstream_control << cur_task_index << " >= " << all_tasks_count << endl;
				sstream_control << "skip" << endl;
				break;
			}
			SendPredictTask( ProcessListNumber, current_status.MPI_SOURCE, cur_task_index, cur_decomp_set_index );
		}
	} // while ( solved_tasks_count < tasks_count )
	
	sstream_control << "after while (solved_tasks_count < all_tasks_count) loop" << endl;
	sstream_control << "solved_tasks_count " << solved_tasks_count << endl;
	sstream_control << "global_deep_point_index " << global_deep_point_index << endl;
	sstream_control << "total_decomp_set_count " << global_deep_point_index << endl;

	if ( verbosity > 2 )
		cout << sstream_control.str() << endl;
	
	GetPredict();
	if ( verbosity > 2 )
		std::cout << "After GetPredict()" << std::endl;

	WritePredictToFile( all_skip_count, whole_time_sec );
	
	if ( verbosity > 2 )
		std::cout << "After WritePredictToFile()" << std::endl;
	
	return true;
}

bool MPI_Predicter :: solverSystemCalling( vec<Lit> &dummy )
{	
	// for every literal make oneliteral clause

	if ( ( rank == 1 ) && ( verbosity > 0 ) )
		std::cout << "start of solverSystemCalling()" << std::endl;
	
	oneliteral_string_vec.resize( dummy.size() );
	std::stringstream sstream;
	for( int i=0; i < dummy.size(); i++ ) {
		sstream << (sign(dummy[i]) ? "-" : "") << var(dummy[i]) + 1;
		oneliteral_string_vec[i] = "\n" + sstream.str() + " 0";
		sstream.str(""); sstream.clear();
	}

	if ( ( rank == 1 ) && ( verbosity > 2 ) )
		std::cout << "oneliteral_string_vec.size() " << oneliteral_string_vec.size() << std::endl;
	
#ifndef WIN32
	truncate( tmp_cnf_process_name.c_str(), template_cnf_size );
#endif
	fstream ofile;
	ofile.open( tmp_cnf_process_name, std::ios::in | std::ios::out );
	// write head of CNF
	sstream << "p cnf " << var_count << " " << clause_count + oneliteral_string_vec.size();
	std::string str = sstream.str();
	unsigned size = str.size();
	if ( size > 20 ) {
		std::cerr << "size > 20" << std::endl;
		return false;
	}
	for ( unsigned i=0; i < 20-size; i++ )
		str += ' '; // add additional wymbols to replace whole 1st line
	ofile.seekp( 0, ofile.beg );
	ofile << str;
	ofile.clear();
	//ofile.close();
	//ofile.open( tmp_cnf_process_name, std::ios::out | std::ios::app );
	ofile.seekp( 0, ofile.end ); // old oneliteral clauses were truncated previosly
	for( auto &x : oneliteral_string_vec )
		ofile << x;
	ofile.close();
	
	string system_str = make_solver_launch_str( solver_name, tmp_cnf_process_name, max_solving_time );
	
	if ( ( rank == 1 ) && ( verbosity > 2 ) )
		std::cout << "system_str " << system_str << std::endl;
	
	fstream current_cnf_out;
	current_cnf_out.open( current_cnf_out_name, std::ios_base :: out );
	//std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	cnf_time_from_node = MPI_Wtime();
	current_cnf_out << Addit_func::exec( system_str );
	//std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	//cnf_time_from_node = time_span.count();
	cnf_time_from_node = MPI_Wtime() - cnf_time_from_node;
	if ( cnf_time_from_node > SOLVER_PARSE_SIMP_TIME )
		cnf_time_from_node -= SOLVER_PARSE_SIMP_TIME;
	if ( cnf_time_from_node <= 0.0 ) {
		std::cerr << "cnf_time_from_node <= 0.0: " << cnf_time_from_node << std::endl;
		return false;
	}
	current_cnf_out.clear(); current_cnf_out.close();
	current_cnf_out.open( current_cnf_out_name.c_str(), std::ios_base :: in );
	
	process_sat_count = 0;
	unsigned str_count = 0;
	while ( getline( current_cnf_out, str ) ) {
		if ( str.find("SATISFIABLE") != std::string::npos )
			process_sat_count = 1;
		str_count++;
	}
	current_cnf_out.close();
	
	if ( !str_count ) {
		std::cerr << "empty solver out file " << current_cnf_out_name << std::endl;
		return false;
	}
	
	return true;
}

bool MPI_Predicter :: solverProgramCalling( vec<Lit> &dummy ) 
{
	int prev_sat_sample_index = -1;
	uint64_t prev_starts, prev_conflicts;
	lbool ret;
	stringstream sstream;
	string cur_state_file_name;
			
	if ( ( verbosity > 2 ) && ( rank == 1 ) )
		std::cout << "Before getting prev stats from Solver" << std::endl;
	
	// make Solver for every problem to make them independent from each other
	Solver *S;
	S = new Solver();
	S->addProblem(cnf);
	S->pdsat_verbosity  = verbosity;
	S->IsPredict        = IsPredict;
	S->max_solving_time = max_solving_time;
	S->rank             = rank;
	S->core_len         = core_len;
	S->start_activity   = start_activity;
	if ( solver_name.find("minigolf") != std::string::npos )
		S->cur_hack_type = hack_minigolf;
	
	prev_starts    = S->starts;
	prev_conflicts = S->conflicts;
	isSolvedOnPreprocessing = 0;
			
	S->resetVarActivity();

	cnf_time_from_node = MPI_Wtime( );
	if ( ( verbosity > 2 ) && ( rank == 1 ) )
		std::cout << "Before S->solveLimited( dummy )" << std::endl;
	ret = S->solveLimited( dummy );
	if ( ( verbosity > 2 ) && ( rank == 1 ) )
		std::cout << "After S->solveLimited( dummy )" << std::endl;
	if ( evaluation_type == "time" )
		cnf_time_from_node = MPI_Wtime( ) - cnf_time_from_node;
	else if ( evaluation_type == "propagation" )
		cnf_time_from_node = (double)S->propagations;
			
	if ( ( S->starts - prev_starts <= 1 ) && ( S->conflicts == prev_conflicts ) )
		isSolvedOnPreprocessing = 1;  // solved by BCP
			
	if ( ( te > 0 ) && ( ret == l_False ) ) { // in ro es te mode all instances are satisfiable
		std::cerr << "( te > 0 ) && ( ret == l_False ) " << std::endl;
		exit(1);
	}
			
	S->getActivity( full_var_choose_order, var_activity, activity_vec_len ); // get activity of Solver

	if ( ( verbosity > 2 ) && ( rank == 1 ) )
		std::cout << "After S->getActivity" << std::endl;
	if ( cnf_time_from_node < MIN_SOLVE_TIME ) // TODO. maybe 0 - but why?!
		cnf_time_from_node = MIN_SOLVE_TIME;
	if ( ret == l_True ) {
		process_sat_count++;
		if ( te == 0 ) {
 			std::cout << "SAT found" << std::endl;
			std::cout << "process_sat_count " << process_sat_count << std::endl;
			b_SAT_set_array.resize( S->model.size() );
			for ( int i=0; i < S->model.size(); i++ )
				b_SAT_set_array[i] = ( S->model[i] == l_True) ? 1 : 0;
			if ( !AnalyzeSATset( ) ) { 	// check res file for SAT set existing
				std::cout << "Error in Analyzer procedute" << std::endl;
				return false;
			}
		}
		/*if ( var_choose_order.size() <= 30 ) { // write activity of hard SAT instance
			S->getActivity( all_vars_set, all_var_activity, all_vars_set.size() );
			ofstream ofile( "bkv_activity" );
			for ( unsigned i=0; i < all_vars_set.size(); i++)
				ofile << all_var_activity[i] << " ";
			ofile.close();
		}
		if ( var_choose_order.size() <= blob_var_count ) {
			sstream << "blob_point_set_len_" << var_choose_order.size() << "_rank_" << rank << "_num_" << cur_point_number;
			cur_state_file_name = sstream.str();
			sstream.clear(); sstream.str("");
			S->clearDB();
			S->saveState( cur_state_file_name );
			sstream << "blob_point_set_len_" << var_choose_order.size() << "_rank_" << rank << "_num_" << cur_point_number << "_set";
			cur_state_file_name = sstream.str();
			sstream.clear(); sstream.str("");
			ofstream ofile( cur_state_file_name.c_str() );
			for ( unsigned i=0; i < var_choose_order.size(); i++ )
				ofile << var_choose_order[i] << " ";
			ofile.close();
		}*/
	}
	delete S;
	/*S->clearDB();
	S->clearPolarity();
	S->clearParams();*/

	return true;
}

bool MPI_Predicter :: ComputeProcessPredict( )
{	
	MPI_Status status;
	int message_size;
	stringstream sstream;
	// get core_len before getting tasks
	MPI_Recv( &var_count,        1, MPI_UNSIGNED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
	MPI_Recv( &core_len,         1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
	MPI_Recv( &activity_vec_len, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
	MPI_Recv( &known_last_bits,  1, MPI_UNSIGNED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
	MPI_Recv( &input_var_num,    1, MPI_UNSIGNED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
	MPI_Recv( &start_activity,   1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
	if ( input_var_num == 0 ) {
		std::cerr << "input_var_num == 0" << endl;
		exit(1);
	}
	if ( rank == 1 ) {
		std::cout << "rank 1" << endl;
		std::cout << "var_count "        << var_count << std::endl;
		std::cout << "core_len "         << core_len << std::endl;
		std::cout << "activity_vec_len " << activity_vec_len << std::endl;
		std::cout << "known_last_bits "  << activity_vec_len << std::endl;
		std::cout << "input_var_num "    << input_var_num << std::endl;
		std::cout << "start_activity "   << start_activity << std::endl;
	}
	int *full_local_decomp_set = new int[core_len];
	MPI_Recv( full_local_decomp_set, core_len, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
	full_var_choose_order.resize( core_len );
	for ( unsigned i=0; i < core_len; ++i )
		full_var_choose_order[i] = full_local_decomp_set[i];
	delete[] full_local_decomp_set;
	all_vars_set.resize( var_count ); // all vars including additional anf keystream
	for ( unsigned i=0; i < all_vars_set.size(); i++ )
		all_vars_set[i] = i+1;
	sstream << "collisions_rank" << rank;
	string collisions_file_name = sstream.str();
	sstream.str(""); sstream.clear();

	if ( isSolverSystemCalling ) {
		// get data and size of template cnf file
		std::fstream file;
		std::string str;
		file.open( input_cnf_name, std::ios_base::in );
		if ( !file.is_open() ) {
			std::cerr << "!file.is_open()" << std::endl;
			return false;
		}
		bool isFirstTemplateStringRead = false;
		clause_count = 0;
		while ( getline( file, str ) ) {
			str.erase( std::remove(str.begin(), str.end(), '\r'), str.end() ); // remove CR symbol
			if ( isFirstTemplateStringRead )
				template_sstream << std::endl;
			if ( str[0] != 'c' ) {
				template_sstream << str;
				isFirstTemplateStringRead = true;
				if ( str[0] != 'p' )
					clause_count++;
			}
		}
		file.close();
		// get final size in bytes
		template_sstream.clear();
		template_sstream.seekg(0, template_sstream.end);
		template_cnf_size = template_sstream.tellg();
		template_sstream.clear();
		
		sstream << rank;
		tmp_cnf_process_name = "./tmp_cnf/tmp_cnf_process_" + sstream.str();
		current_cnf_out_name = "./tmp_cnf/out_solver_" + sstream.str();
		sstream.clear(); sstream.str("");

		// add template file to tmp cnf file of every computing process
		ofstream tmp_cnf_process;
		tmp_cnf_process.open( tmp_cnf_process_name, std::ios_base::out );
		tmp_cnf_process << template_sstream.str();
		tmp_cnf_process.close();
	}
	
	// read file with CNF once
	if ( !isSolverSystemCalling ) {
		minisat22_wrapper m22_wrapper;
		ifstream in( input_cnf_name );
		m22_wrapper.parse_DIMACS_to_problem(in, cnf);
		in.close();
		/*S = new Solver();
		S->addProblem(cnf);
		S->pdsat_verbosity  = verbosity;
		S->IsPredict        = IsPredict;
		S->max_solving_time = max_solving_time;
		S->rank             = rank;
		S->core_len         = core_len;
		S->start_activity   = start_activity;
		if ( solver_name.find("minigolf") != std::string::npos )
			S->cur_hack_type = hack_minigolf;*/
	}
	
	if ( te > 0 ) { // ro es te mode
		MPI_Recv( &first_stream_var_index,  1, MPI_UNSIGNED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		if ( rank == 1 )
			cout << "received first_stream_var_index " << first_stream_var_index << endl;
		int stream_char_len, state_char_len;
		unsigned stream_vec_len, state_vec_len;
		MPI_Probe( 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		MPI_Get_count( &status, MPI_CHAR, &stream_char_len );
		//if ( rank == 1 )
		//	cout << "stream_char_len " << stream_char_len << endl;
		char *stream_arr = new char[stream_char_len];
		MPI_Recv( stream_arr, stream_char_len, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		stream_vec_len = stream_char_len / cnf_in_set_count;
		stream_vec_vec.resize( cnf_in_set_count );
		for( unsigned i=0; i < cnf_in_set_count; i++ )
			for( unsigned j=0; j < stream_vec_len; j++ )
				stream_vec_vec[i].push_back( stream_arr[i*stream_vec_len + j] ? true : false );
		delete[] stream_arr;
		if ( rank == 1 ) {
			cout << "stream_vec_len " << stream_vec_len << endl; 
			cout << "stream_vec_vec.size() " << stream_vec_vec.size() << endl;
		}
		MPI_Probe( 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		MPI_Get_count( &status, MPI_CHAR, &state_char_len );
		if ( rank == 1 )
			cout << "state_char_len " << state_char_len << endl;
		char *state_arr = new char[state_char_len];
		MPI_Recv( state_arr, state_char_len, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		state_vec_len = state_char_len / cnf_in_set_count;
		state_vec_vec.resize( cnf_in_set_count );
		for( unsigned i=0; i < cnf_in_set_count; i++ )
			for( unsigned j=0; j < state_vec_len; j++ )
				state_vec_vec[i].push_back( state_arr[i*state_vec_len + j] ? true : false );
		delete[] state_arr;
		if ( rank == 1 ) {
			cout << "state_vec_len " << state_vec_len << endl;
			cout << "state_vec_vec.size() " << state_vec_vec.size() << endl;
		}
	}
	
	var_choose_order = full_var_choose_order;
	
	if ( rank == 1 ) {
		cout << "Received core_len "         << core_len         << endl;
		cout << "Received activity_vec_len " << activity_vec_len << endl;
		cout << "Received known_last_bits "  << known_last_bits  << endl;
		cout << "Received full_var_choose_order " << endl;
		for ( unsigned i=0; i < full_var_choose_order.size(); ++i )
			cout << full_var_choose_order[i] << " ";
		cout << endl;
	}
	
	current_task_index = 0;
	cnf_time_from_node = 0.0;
	int ProcessListNumber;
	var_activity = new double[activity_vec_len];
	all_var_activity = new double[all_vars_set.size()];
	bool isFirstDecompSetReceived = false, isNewDecompSetReceived = false;
	unsigned large_message_count = 0, small_message_count = 0;
	string polarity_file_name;
	isSolvedOnPreprocessing = 0;
	
	for (;;) {
		do {// get index of current task missing stop-messages
			MPI_Probe( 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
			MPI_Get_count( &status, MPI_INT, &message_size );
			if ( message_size == 0 ) {
				cerr << "After MPI_Get_count( &status, MPI_INT, &message_size ) message_size " << message_size << endl;
				return false;
			}
			if ( ( verbosity > 0 ) && ( rank == 1 ) )
				cout << "recv message_size " << message_size << endl;
			if ( message_size == 1 ) {
				MPI_Recv( &current_task_index, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
				isNewDecompSetReceived = false;
			}
			else if ( message_size > 1 )
				break; // stop and recv array message
			else if ( ( message_size > 1 ) && ( current_task_index == -1 ) ) {
				cerr << "( message_size > 1 ) && ( current_task_index == -1 )" << endl;
				cerr << "message_size " << message_size << endl;
				return false;
			}
			if ( ( current_task_index == -1 ) && ( verbosity > 1 ) && ( rank == -1 ) )
				cout << "on node " << rank << " stop-message was skipped" << endl;
		} while ( current_task_index == -1 ); // skip stop-messages
		
		ProcessListNumber = status.MPI_TAG;
		
		if ( message_size-1 > (int)full_var_choose_order.size() ) {
			cerr << "message_size > full_var_choose_order.size()" << endl;
			cerr << message_size << " > " << full_var_choose_order.size() << endl;
			return false;
		}
		
		if ( ( verbosity > 0 ) && ( rank == 1 ) ) {
			cout << "small_message_count " << small_message_count << endl;
			cout << "large_message_count " << large_message_count << endl;
		}
		
		if ( message_size == 1 )
			small_message_count++;
		else if ( message_size > 1 ) {
			large_message_count++;
			if ( isFirstDecompSetReceived )
				delete[] array_message;
			array_message = new int[message_size];
			MPI_Recv( array_message, message_size, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
			cur_point_number++;
			isNewDecompSetReceived = true;
			if ( ( verbosity > 2 ) && ( rank == 1 ) ) {
				std::cout << "Received array_message" << std::endl;
				for ( int i = 0; i < message_size; ++i )
					cout << array_message[i] << " ";
				std::cout << std::endl;
			}
			current_task_index = array_message[0];
			var_choose_order.resize( message_size - 1 );
			for ( int i=1; i < message_size; ++i )
				var_choose_order[i-1] = array_message[i];
			if ( ( verbosity > 2 ) && ( rank == 1 ) ) {
				std::cout << "Received var_choose_order" << endl;
				for ( unsigned i = 0; i < var_choose_order.size(); ++i )
					std::cout << var_choose_order[i] << " ";
				std::cout << std::endl;
			}
			/*if ( solver_type == 4 ) {
				if ( isFirstDecompSetReceived ) // if not first time, delete old data
					delete S;
			}*/
			isFirstDecompSetReceived = true;
		}
		
		if ( ( verbosity > 0 ) && ( rank == 1 ) )
			cout << "current_task_index " << current_task_index << endl;

		if ( !isFirstDecompSetReceived ) {
			cerr << "!IsFirstDecompSetReceived before solving" << endl;
			return false;
		}

		// get array of assignments
		vec<Lit> dummy;
		unsigned cur_var_ind;
		int sat_sample_index;
		int cur_stream_index;
		int val;
		
		if ( te > 0 ) { // make satisfiable instanse by known state and keystream
			dummy.clear();
			sat_sample_index = current_task_index % cnf_in_set_count;
			/*if ( sat_sample_index == prev_sat_sample_index ) {
				S->clearPolarity(); // if same keystrem, clear phase saving
				sstream << rank;
				polarity_file_name = "polarity_file_name_rank" + sstream.str();
				sstream.clear(); sstream.str("");
				ofstream ofile( polarity_file_name, ios_base :: out | ios_base :: app );
				ofile << "sat_sample_index == prev_sat_sample_index" << endl;
				ofile << sat_sample_index << " == " << prev_sat_sample_index << endl;
				ofile << "prev_set_size " << prev_set_size << endl;
				ofile << "cur_set_size " << var_choose_order.size() << endl;
				ofile << endl;
				ofile.close();
			}
			prev_sat_sample_index = sat_sample_index;
			prev_set_size = var_choose_order.size();*/
				
			for ( auto &x : var_choose_order ) {
				cur_var_ind = x-1;
				if ( cur_var_ind > state_vec_vec[sat_sample_index].size()-1 ) {
					cerr << "cur_var_ind > state_vec_vec[sat_sample_index].size()-1" << endl;
					exit(1);
				}
				dummy.push( (state_vec_vec[sat_sample_index][cur_var_ind]) ? mkLit( cur_var_ind ) : ~mkLit( cur_var_ind ) );
			}
			if ( known_last_bits ) { // add some last known bits
				for ( unsigned i=0; i < known_last_bits; i++ ) {
					cur_var_ind = core_len - known_last_bits + i;
					dummy.push( (state_vec_vec[sat_sample_index][cur_var_ind]) ? mkLit( cur_var_ind ) : ~mkLit( cur_var_ind ) );
				}
			}
			cur_stream_index = 0;
			for ( auto it = stream_vec_vec[sat_sample_index].begin(); it != stream_vec_vec[sat_sample_index].end(); it++ ) {
				cur_var_ind = first_stream_var_index + cur_stream_index;
					dummy.push( *it ? mkLit( cur_var_ind ) : ~mkLit( cur_var_ind ) );
				cur_stream_index++;
			}
		}
		else { // make random values of decomp variables
			dummy.resize( var_choose_order.size() );
			for ( unsigned i=0; i < var_choose_order.size(); ++i ) {
				val = var_choose_order[i] - 1;
				dummy[i] = bool_rand( gen ) ? mkLit( val ) : ~mkLit( val );
			}
		}
			
		if ( ( verbosity > 2 ) && ( rank == 1 ) ) {
			std::cout << std::endl;
			std::cout << "dummy size " << dummy.size() << std::endl;
			for ( int i=0; i < dummy.size(); ++i )
				std::cout << dummy[i].x << " ";
			std::cout << std::endl;
		}
		
		cnf_time_from_node = 0.0;
		process_sat_count = 0;
		if ( !isSolverSystemCalling ) {
			if ( !solverProgramCalling( dummy ) )
				MPI_Abort( MPI_COMM_WORLD, 0 );
		}
		else {
			if ( !solverSystemCalling( dummy ) )
				MPI_Abort( MPI_COMM_WORLD, 0 );
		}
		
		MPI_Send( &current_task_index,         1, MPI_INT,    0, ProcessListNumber, MPI_COMM_WORLD );
		MPI_Send( &process_sat_count,          1, MPI_INT,    0, ProcessListNumber, MPI_COMM_WORLD );
		MPI_Send( &cnf_time_from_node,         1, MPI_DOUBLE, 0, ProcessListNumber, MPI_COMM_WORLD );
		if ( !isSolverSystemCalling ) {
			MPI_Send( &isSolvedOnPreprocessing,    1, MPI_INT,    0, ProcessListNumber, MPI_COMM_WORLD );
			MPI_Send( var_activity, activity_vec_len, MPI_DOUBLE, 0, ProcessListNumber, MPI_COMM_WORLD );
		}
		
		if ( verbosity > 0 )
			std::cout << "rank " << rank << " sended decision" << std::endl;
	}
	
	delete[] var_activity;
	delete[] all_var_activity;
 	delete[] array_message;
	MPI_Finalize();
	return true;
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
	
	if ( known_point_file.is_open() ) { // get known point
		cout << "known_point_file opened " << known_point_file_name << endl;
		int ival;
		var_choose_order.resize(0);
		while ( known_point_file >> ival )
			var_choose_order.push_back( ival );
		best_var_num = var_choose_order.size(); 
		known_point_file.close();
	}
	else {
		sstream << "schema_type " << schema_type << endl;
		if ( ( schema_type != "rand" ) && ( var_choose_order.size() == 0 ) ) { // if schema_type was set by user
			full_mask_var_count = predict_to;
			MakeVarChoose();
		}
		else if ( schema_type == "rand" ) { // if no file with known point then get random init point
			vector<unsigned> rand_arr;
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
	}
	sort( var_choose_order.begin(), var_choose_order.end() );
	sstream << "var_choose_order" << endl;
	for ( unsigned i = 0; i < var_choose_order.size(); i++ )
		sstream << var_choose_order[i] << " ";

	unsigned array_message_size = var_choose_order.size() + 1;
	array_message = new int[array_message_size];
	for ( unsigned i=1; i < array_message_size; ++i )
		array_message[i] = var_choose_order[i-1];
	
	sstream << endl << endl;
	deep_predict_file.open( deep_predict_file_name.c_str(), ios_base::out | ios_base::app );
	deep_predict_file << sstream.rdbuf();
	deep_predict_file.close();
}

bool MPI_Predicter :: DeepPredictFindNewUncheckedArea( stringstream &sstream ) 
{
	string str;
	unsigned L1_erased_count = 0;
	unsigned L2_erased_count = 0;
	list<unchecked_area> :: iterator L2_it;
	bool isChoosingByActivity;
	
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
		boost::dynamic_bitset<> bs = IntVecToBitsetPredict( var_choose_order );
		bool IsCorrectAddingArea = false;
		// remove points from L2 that too far from record
		// and find L2 for current record
		L2_it = L2.begin();
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
		str = string( str.rbegin(), str.rend() );
		sstream << "current_unchecked_area center " << endl << str << endl;
		to_string( current_unchecked_area.checked_points, str );
		str = string( str.rbegin(), str.rend() );
		sstream << "current_unchecked_area checked_points " << endl << str << endl;;
		// make initial values
		cur_vars_changing = 1; // start again from Hamming distance == 1
		unupdated_count = 0;
	} else { // if there were no better points in checked area
		sstream << endl << "---Record not updated---" << endl;
		checked_area c_a;
		unupdated_count++;
		sstream << "unupdated_count " << unupdated_count << endl;
		if ( deep_predict >= 5 ) { // tabu search mode
			boost::dynamic_bitset<> bs, xor_bs;
			unsigned min_hamming_distance;
			bool IsAdding;
			unsigned rand_power_value, rand_L2_start_search_index;
			list<unchecked_area> L2_matches;
			vector<point_struct> point_struct_vec;
			vector<int> power_values;
			unsigned L2_index;
			bs = IntVecToBitsetPredict( real_var_choose_order ); // current center of searching
			unsigned max_checked = 0;
			// find new unchecked area
			sstream << "finding new unchecked_area" << std::endl;
			sstream << "ts_strategy " << ts_strategy << std::endl;
			if ( verbosity > 0 )
				std::cout << "finding new unchecked_area" << std::endl;
			// find min hamming distance of points in L2 from current record point
			min_hamming_distance = (unsigned)core_len;
			for ( L2_it = L2.begin(); L2_it != L2.end(); ++L2_it ) {
				xor_bs = (*L2_it).center ^ bs;
				if ( xor_bs.count() == 0 ) {
					std::cerr << "xor_bs == 0. current center in L2" << std::endl;
					std::cerr << "(*L2_it).checked_points.count() " << (*L2_it).checked_points.count() << std::endl;
					MPI_Abort( MPI_COMM_WORLD, 0 );
				}
				if ( xor_bs.count() < min_hamming_distance )
					min_hamming_distance = xor_bs.count();
			}
			sstream << "min hamming distance from L2 " << min_hamming_distance << std::endl;
			if ( min_hamming_distance > max_L2_hamming_distance ) {
				std::cerr << "min_hamming_distance > max_L2_hamming_distance " << std::endl;
				std::cerr << min_hamming_distance << " > " << max_L2_hamming_distance << std::endl;
				return false;
			}
			/*if ( cur_vars_changing != min_hamming_distance ) {
				cur_vars_changing = min_hamming_distance;
				sstream << "cur_vars_changing changed to " << cur_vars_changing << endl;
			}*/
			// remember matches
			for ( L2_it = L2.begin(); L2_it != L2.end(); L2_it++ ) {
				xor_bs = (*L2_it).center ^ bs;
				if ( xor_bs.count() == min_hamming_distance )
					L2_matches.push_back( *L2_it );
			}
			
			if ( L2_matches.size() == 0 ) {
				cerr << "Error. L2_matches.size() == 0" << endl; exit(1);
			}
			if ( verbosity > 0 )
				cout << "L2_matches.size() " << L2_matches.size() << endl;
			sstream << "L2_matches.size() " << L2_matches.size() << endl;
			isChoosingByActivity = true;
			
			if ( isChoosingByActivity ) {
				sstream << "isChoosingByActivity " << isChoosingByActivity << endl;
				switch ( ts_strategy ) { // randomly choose weight of point from L2 and go to such random point
					case 0:
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
						rand_L2_start_search_index = uint_rand( gen ) % L2_matches.size();
						rand_power_value           = uint_rand( gen ) % power_values.size();
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
						// try to find to another side if we havn't found point
						if ( !IsAdding ) {
							L2_index = 0;
							for ( L2_it = L2_matches.begin(); L2_it != L2_matches.end(); ++L2_it ) {
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
						break;
					case 1: // sort areas from L2 by  median of var activities of centers
						// update activity of points from L2
						for ( L2_it = L2_matches.begin(); L2_it != L2_matches.end(); ++L2_it ) {
							(*L2_it).med_var_activity = 0;
							for ( unsigned i=0; i < (*L2_it).center.size(); ++i )
								if ( (*L2_it).center[i] )
									(*L2_it).med_var_activity += total_var_activity[i];
							(*L2_it).med_var_activity /= (*L2_it).center.count();
						}
						L2_matches.sort( ua_compareByActivity );
						if ( verbosity > 0 ) { 
							cout << "L2 after sorting. total_var_activity : center.count()";
							for ( L2_it = L2_matches.begin(); L2_it != L2_matches.end(); ++L2_it )
								cout << (*L2_it).med_var_activity << " : " << (*L2_it).center.count() << endl;
						}
						current_unchecked_area = (*L2_matches.begin());
						break;
				}
			}
			
			L2_matches.clear();
			
			if ( verbosity > 0 )
				cout << "bofore BitsetToIntVecPredict()" << endl;

			var_choose_order = BitsetToIntVecPredict( current_unchecked_area.center );
			to_string( current_unchecked_area.center, str );
			str = string( str.rbegin(), str.rend() );
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

	map<unsigned, unsigned> point_map; // set count of points + count of such points
	map<unsigned, unsigned> :: iterator point_map_it;
	for ( L2_it = L2.begin(); L2_it != L2.end(); ++L2_it ){
		point_map_it = point_map.find( (*L2_it).center.count() );
		if ( point_map_it == point_map.end() )
			point_map.insert( pair<unsigned,unsigned>( (*L2_it).center.count(), 1) );
		else
			(*point_map_it).second++;
	}
	ofstream ofile("L2_list", ios_base :: out );
	ofile << "L2 point # : set_count : size" << endl;
	map<unsigned, unsigned> :: reverse_iterator point_map_rit;
	unsigned k = 1;
	for ( point_map_rit = point_map.rbegin(); point_map_rit != point_map.rend(); ++point_map_rit )
		ofile << k++ << " : " << (*point_map_rit).first << " : " << (*point_map_rit).second << endl;
	ofile.close();
	point_map.clear();
	
	return true;
}

bool MPI_Predicter :: DeepPredictMain( )
{
	// make renadom init point for 1st iteration
	// for other ones as init use best point of prev iteration
	int ProcessListNumber = 0;
	int stop_message = -1;
	fstream deep_predict_file;
	stringstream sstream;
	whole_deep_time = Minisat::cpuTime();
	stringstream sstream_control;

	cout << "DeepPredictMain() started" << endl;

	// for equal conditions for every dimension
	global_count_var_changing.resize( max_var_deep_predict );
	for ( unsigned i = 0; i < global_count_var_changing.size(); i++ )
		global_count_var_changing[i] = 0;
	
	// read from file if one exists. if not, get random vector
	GetInitPoint( );
	double current_time = 0;
	string str;
	
	while (
		   ( ( deep_predict == 5 ) && 
		     ( ( cur_temperature == 0 ) ||
		       ( cur_temperature > min_temperature ) ) 
			 )
			||
			( ( deep_predict == 3 ) || ( deep_predict == 4 ) || ( deep_predict >= 6 ) )  
		   )	
	{
		if ( isFirstPoint )
		    current_predict_start_time = MPI_Wtime( );
		
		current_time = MPI_Wtime( );
		if ( !GetDeepPredictTasks() )
		{ cerr << "Error in GetDeepPredictTasks" << endl; return false; }
		if ( verbosity > 0 )
			cout << " GetDeepPredictTasks() done" << endl;
		whole_get_deep_tasks_time += MPI_Wtime( ) - current_time;
		//cout << "GetDeepPredictTasks() time " << MPI_Wtime( ) - current_time << endl;

		current_time = MPI_Wtime( );
		if ( !PrepareForPredict( ) )
			{ cerr << "Error in PrepareForPredict" << endl; return false; }
		if ( verbosity > 0 )
			cout << "PrepareForPredict() done" << endl;

		if ( !ControlProcessPredict( ProcessListNumber++, sstream_control ) ) {
			MPI_Abort( MPI_COMM_WORLD, 0 );
		}

		if ( verbosity > 0 )
			cout << "ControlProcessPredict() done" << endl;

		deep_predict_file.open( deep_predict_file_name.c_str(), ios_base::out | ios_base::app );
		sstream.str( "" ); sstream.clear( );
		
		sstream << sstream_control.str();
		
		if ( ( isFirstPoint ) && ( global_deep_point_index == total_decomp_set_count ) ) {
			sstream << "First point" << endl;
			// set new unchecked area
			current_unchecked_area = *L2.begin();
			to_string( current_unchecked_area.center, str );
			str = string( str.rbegin(), str.rend() );
			sstream << "current_unchecked_area center " << endl << str << endl;
			to_string( current_unchecked_area.checked_points, str );
			str = string( str.rbegin(), str.rend() );
			sstream << "current_unchecked_area checked_points " << endl << str << endl;
			isFirstPoint = false;
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
		
		// stop all current tasks if not first stage (in this stage all problems must be solved)
		/*if ( !IsFirstStage ) {
			if ( verbosity > 1 )
				cout << "Extra stop sending " << endl;
			for ( int i = 1; i < corecount; i++ ) {
				MPI_Isend( &stop_message, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &request );
				if ( verbosity > 0 )
					cout << "stop-message was send to node # " << i << endl;
			}
		}*/
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
	sstream << "***whole_deep_time "           << Minisat::cpuTime() - whole_deep_time << " s" << endl;
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
	vector<int> ::iterator vec_it;
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &corecount );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	
	if ( ( solver_name.find( "." ) != std::string::npos ) ||
		 ( solver_name.find( "/" ) != std::string::npos )
		 )
	{
		isSolverSystemCalling = true;
		ts_strategy = 0;
	}
	
	IsPredict = true;
	if ( corecount < 2 ) { 
		std::cerr << "Error. corecount < 2" << std::endl; return false; 
	}
	
	vec_IsDecompSetSendedToProcess.resize( corecount - 1 );
	if ( te > 0 )
		max_solving_time = te; // for (ro, es, te) predict mode

	int stream_char_len, state_char_len;
	
	if ( rank == 0 ) { // control node
		std::cout << "MPI_Predict is running " << std::endl;
		
		if ( !ReadIntCNF( ) ) { 
			std::cerr << "Error in ReadIntCNF" << std::endl; return 1; 
		}
		
		if ( predict_to > MAX_CORE_LEN ) {
			std::cout << "Warning. predict_to > MAX_PREDICT_TO. Changed to MAX_PREDICT_TO" << std::endl;
			std::cout << predict_to << " > " << MAX_CORE_LEN << std::endl;
			predict_to = MAX_CORE_LEN;
		}
		if  ( predict_to > (int)core_len ) {
			std::cout << "core_len changed to predict_to" << std::endl;
			std::cout << core_len << " changed to " << predict_to << std::endl;
			core_len = predict_to;
		}

		if ( ( !predict_to ) && ( core_len ) ) {
			std::cout << "predict_to changed to core_len" << std::endl;
			std::cout << predict_to << " changed to " << core_len << std::endl;
			predict_to = core_len;
		}

		activity_vec_len = core_len;
		var_activity = new double[activity_vec_len];
		for ( unsigned i=0; i < activity_vec_len; i++ )
			var_activity[i] = 0.0;
		total_var_activity.resize( activity_vec_len );
		//for( auto &x : total_var_activity ) x = 0;
		for( vector<double> :: iterator it = total_var_activity.begin(); it != total_var_activity.end(); ++it )
			*it = 0;
		
		int *full_local_decomp_set = new int[core_len];
		for( unsigned i=0; i < core_len; ++i )
			full_local_decomp_set[i] = full_var_choose_order[i];
		
		// send core_len once to every compute process
		for( int i=0; i < corecount-1; ++i ) {
			MPI_Send( &var_count,        1, MPI_UNSIGNED, i + 1, 0, MPI_COMM_WORLD );
			MPI_Send( &core_len,         1, MPI_INT,  i + 1, 0, MPI_COMM_WORLD );
			MPI_Send( &activity_vec_len, 1, MPI_INT,  i + 1, 0, MPI_COMM_WORLD );
			MPI_Send( &known_last_bits,  1, MPI_UNSIGNED, i + 1, 0, MPI_COMM_WORLD );
			MPI_Send( &input_var_num,    1, MPI_UNSIGNED, i + 1, 0, MPI_COMM_WORLD );
			MPI_Send( &start_activity,   1, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD );
			MPI_Send( full_local_decomp_set, core_len, MPI_INT,  i + 1, 0, MPI_COMM_WORLD );
		}
		delete[] full_local_decomp_set;
		
		if ( te > 0 ) {
			first_stream_var_index = var_count - keystream_len;
			cout << "first_stream_var_index " << first_stream_var_index << endl;

			MakeSatSample( state_vec_vec, stream_vec_vec );
			stream_char_len = stream_vec_vec.size() * stream_vec_vec[0].size();
			state_char_len  = state_vec_vec.size()  * state_vec_vec[0].size();
			cout << "stream_vec_vec.size() " << stream_vec_vec.size() << endl;
			cout << "state_vec_vec.size() "  << state_vec_vec.size()  << endl;
			cout << "stream_char_len "       << stream_char_len       << endl;
			cout << "state_char_len "        << state_char_len        << endl;
			char *stream_arr = new char[stream_char_len];
			char *state_arr  = new char[state_char_len];
			int k = 0;
			for ( auto &x : stream_vec_vec )
				for ( auto y = x.begin(); y != x.end(); y++ )
					stream_arr[k++] = *y ? 1 : 0;
			cout << "stream_arr size " << k << endl;
			if ( k != stream_char_len ) {
				cerr << "k != stream_char_len" << endl;
				cerr << k << " != " << stream_char_len << endl;
				exit(1);
			}
			k = 0;
			for ( auto &x : state_vec_vec )
				for ( auto y = x.begin(); y != x.end(); y++ )
					state_arr[k++] = *y ? 1 : 0;
			cout << "state_arr size " << k << endl;
			if ( k != state_char_len ) {
				cerr << "k != state_char_len" << endl;
				cerr << k << " != " << state_char_len << endl;
				exit(1);
			}
			
			for( int i=0; i < corecount-1; ++i ) {
				MPI_Send( &first_stream_var_index,  1, MPI_UNSIGNED, i + 1, 0, MPI_COMM_WORLD );
				MPI_Send( stream_arr, stream_char_len, MPI_CHAR,     i + 1, 0, MPI_COMM_WORLD );
				MPI_Send( state_arr,  state_char_len,  MPI_CHAR,     i + 1, 0, MPI_COMM_WORLD );
			}
			cout << "stream_arr and state_arr sended" << endl;
			
			delete[] stream_arr;
			delete[] state_arr;
		}
		
		array_message = NULL;
		
		std::cout << "verbosity "     << verbosity           << std::endl;
		std::cout << "solver_name "   << solver_name         << std::endl;
		cout << "schema_type "   << schema_type         << endl;
		cout << "corecount "     << corecount           << endl;
		cout << "deep_predict  " << deep_predict        << endl;
		cout << "predict_from "  << predict_from        << endl;
		cout << "predict_to "    << predict_to          << endl;
		cout << "cnf_in_set_count " << cnf_in_set_count << endl; 
		cout << "proc_count "    << proc_count          << endl;
		cout << "core_len "      << core_len            << endl;
		cout << "input_var_num " << input_var_num       << endl;
 		cout << "start_activity "    << start_activity << endl;
		cout << "max_var_deep_predict " << max_var_deep_predict << endl;
		cout << "start_temperature_koef " << start_temperature_koef << endl;
		cout << "point_admission_koef " << point_admission_koef << endl;
		cout << "ts_strategy " << ts_strategy << endl;
		cout << "IsFirstStage " << IsFirstStage << endl;
		cout << "max_L2_hamming_distance " << max_L2_hamming_distance << endl;
		cout << "evaluation_type " << evaluation_type << endl;
		cout << "te for (ro, es, te) " << te << endl;
		cout << "er for (ro, es, te) " << er << endl;
		//cout << "penalty for (ro, es, te) " << penalty << endl;
		cout << "keystream_len " << keystream_len << endl;
		std::cout << "blob_var_count " << blob_var_count << endl;
		std::cout << "isSolverSystemCalling " << isSolverSystemCalling << std::endl;
		std::cout << std::endl;
		
		DeepPredictMain();
		
		delete[] var_activity;
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
bool MPI_Predicter :: PrepareForPredict()
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
	cnf_issat_arr.resize( all_tasks_count );
	cnf_prepr_arr.resize( all_tasks_count );
	
	for( unsigned i = 0; i < all_tasks_count; i++ ) {
		cnf_status_arr[i]     = 0; // status WAIT
		cnf_start_time_arr[i] = 0.0;
		cnf_real_time_arr[i]  = 0.0;
		cnf_prepr_arr[i]      = 0;
		cnf_issat_arr[i]      = false;
	}

	cout << "all_tasks_count " << all_tasks_count << endl;
	
	// count of solved subproblems from sample with time limit
	solved_in_time_arr.resize( decomp_set_arr.size() );
	// 
	time_limit_arr.resize( decomp_set_arr.size() );
	// sum times
	sum_time_arr.resize( decomp_set_arr.size() );
	// array of median times of CNF in set
	med_time_arr.resize( decomp_set_arr.size() );
	// array of predict times
	predict_time_arr.resize( decomp_set_arr.size() );
	// array of partly predict times
	predict_part_time_arr.resize( decomp_set_arr.size() );
	// array of CNF set status
	set_status_arr.resize( decomp_set_arr.size() );
	stopped_cnf_count_arr.resize( decomp_set_arr.size() );
	skipped_cnf_count_arr.resize( decomp_set_arr.size() );
	solved_cnf_count_arr.resize( decomp_set_arr.size() );

	// array of indexes of array of lengths of sets for modelling
	set_index_arr.resize( decomp_set_arr.size() + 1 );
	set_index_arr[0] = 0;
	val = 0;
	// set_index_arr = array of CNF set indexes, it depends on set_len_arr
	for( unsigned i = 0; i < decomp_set_arr.size(); ++i ) {
		val += set_len_arr[i];
		set_index_arr[i + 1] = val;
	}
	
	for( unsigned i = 0; i < decomp_set_arr.size(); ++i ) {
		set_status_arr[i]        = 0;
		stopped_cnf_count_arr[i] = 0;
		skipped_cnf_count_arr[i] = 0;
		solved_cnf_count_arr[i]  = 0;
		sum_time_arr[i]          = 0.0;
		med_time_arr[i]          = 0.0;
		predict_time_arr[i]      = 0.0;
		predict_part_time_arr[i] = 0.0;
		solved_in_time_arr[i]    = 0;
		time_limit_arr[i]        = 0.0;
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
	best_var_num = decomp_set_arr[set_index].var_choose_order.size();
	var_choose_order = decomp_set_arr[set_index].var_choose_order;
	real_var_choose_order = var_choose_order;
	sort( var_choose_order.begin(), var_choose_order.end() );
	
	double last_predict_record_time = MPI_Wtime() - current_predict_start_time;
	current_predict_time += last_predict_record_time;
	current_predict_start_time = MPI_Wtime(); // update time
	
	sstream << "best_predict_time " << best_predict_time << " s" << endl;

	if ( deep_predict == 5 ) { // Simulated annealing
		if ( isSimulatedGranted ) {
			sstream << "*** Simulated granted" << endl;
			//sstream << "rand_num " << rand_num << endl;
			sstream << "delta " << delta << endl;
			sstream << "exp_value " << exp_value << endl;
			isSimulatedGranted = false;
		}
		// init or update temperature if == 0 or > limit
		if ( ( cur_temperature == 0.0 ) || 
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

	sstream << "best_var_num " << best_var_num << std::endl;
	sstream << "best_var_choose_order " << std::endl;
	for ( int j = 0; j < best_var_num; j++ )
		sstream << var_choose_order[j] << " ";
	sstream << std::endl;
	sstream << "global_count_var_changing" << endl;
	for ( int i = 0; i < max_var_deep_predict; i++ )
		sstream << i + 1 << ":" << global_count_var_changing[i] << " ";
	sstream << std::endl << std::endl;

	if ( IsDecDecomp ) {
		sstream << "real best predict time " << real_best_predict_time << " s" << std::endl;
		sstream << "real best var num " << real_best_var_num << std::endl;
		sstream << "real var choose order " << std::endl;
		for ( int j = 0; j < real_best_var_num; j++ )
			sstream << real_var_choose_order[j] << " ";
		sstream << std::endl;
	}
	
	/*if ( 
		( ( deep_predict == 3 ) || ( deep_predict == 5 ) || ( deep_predict == 6 ) )
		&& ( !IsDecDecomp ) && ( !isFirstPoint ) && ( !IsFirstStage )
	   )
		IsRestartNeeded = true;*/

	deep_predict_file << sstream.rdbuf();
	deep_predict_file.close();

	sstream.clear(); sstream.str("");
	record_count++;
	sstream << "predict_" << record_count;
	predict_file_name = sstream.str();
	
	//if ( ( !IsFirstStage ) || ( record_count == 1 ) ) // don't write in first stage - it's expensive
	WritePredictToFile( 0, 0 );
	predict_file_name = "predict";
	
	ofstream graph_file, var_activity_file;
	if ( isFirstPoint ) {
		graph_file.open( "graph_file", ios_base :: out ); // erase info from previous launches 
		if ( te == 0.0 )
			graph_file << "# best_var_num best_predict_time best_sum_time cnf_in_set_count last_predict_record_time current_predict_time";
		else
			graph_file << "# best_var_num best_predict_time solved_in_time cnf_in_set_count last_predict_record_time current_predict_time";
		if ( deep_predict == 5 ) // simulated annealing
			graph_file << " cur_temperature";
		graph_file << std::endl;
		if ( !isSolverSystemCalling ) 
			var_activity_file.open( var_activity_file_name.c_str(), ios_base :: out );
	}
	else {
		graph_file.open( "graph_file", ios_base :: app );
		if ( !isSolverSystemCalling )
			var_activity_file.open( var_activity_file_name.c_str(), ios_base :: app );
	}

	if ( !isSolverSystemCalling )
		var_activity_file << record_count << endl;
	//for( auto &x : total_var_activity )
	//	var_activity_file << x << " ";
	if ( !isSolverSystemCalling ) {
		for( vector<double> :: iterator it = total_var_activity.begin(); it != total_var_activity.end(); ++it )
			var_activity_file << *it << " ";
		var_activity_file << std::endl;
	}
	
	graph_file << record_count << " " << best_var_num << " " << best_predict_time << " ";
	if ( te == 0.0 )
		graph_file << best_sum_time;
	else
		graph_file << best_solved_in_time;
	graph_file << " " << best_cnf_in_set_count << " " << last_predict_record_time << " " << current_predict_time;
	if ( deep_predict == 5 ) 
		graph_file << " " << cur_temperature;

	if ( ( IsFirstStage ) && ( !isFirstPoint ) && ( best_var_num > old_best_var_num ) ) {
		IsFirstStage = false;
		std::cout << "IsFirstStage "      << IsFirstStage      << std::endl;
		std::cout << "best_var_num "      << best_var_num      << std::endl;
		std::cout << "best_predict_time " << best_predict_time << std::endl;
		sstream << std::endl << " *** First stage done ***" << best_predict_time << std::endl;
		graph_file << " first stage done";
	}
	graph_file << " time_limit=" << best_time_limit << std::endl;
	
	graph_file.close();
	if ( !isSolverSystemCalling )
		var_activity_file.close();
}

bool MPI_Predicter :: checkSimulatedGranted( double predict_time )
{
	if ( deep_predict != 5 )
		return false;
	delta = predict_time - best_predict_time;
	if ( delta < 0 )
		return true;
	double rand_num = 0; // get value from [0.1..0.9] with step 0.1
	while ( rand_num == 0 ) {
		rand_num = ( unsigned )( uint_rand( gen ) % 10 );
		rand_num *= 0.1;
	}
	exp_value = exp(-delta / cur_temperature);
	if ( rand_num < exp_value ) {
		isSimulatedGranted = true;	
		return true;
	}
	else return false;
}

vector<int> MPI_Predicter :: BitsetToIntVecPredict( boost::dynamic_bitset<> &bs )
{
	vector<int> vec_int;
	for ( unsigned i=0; i < bs.size(); ++i )
		if ( bs[i] ) vec_int.push_back( full_var_choose_order[i] );
	return vec_int;
}

boost::dynamic_bitset<> MPI_Predicter :: IntVecToBitsetPredict( vector<int> &variables_vec )
{
	boost::dynamic_bitset<> bs( core_len );
	for ( unsigned i=0; i<variables_vec.size(); ++i )
		bs.set( core_var_indexes.find( variables_vec[i] )->second ); // set element to 1
	return bs;
}

double MPI_Predicter :: getCurPredictTime( unsigned cur_var_num, int cur_cnf_in_set_count, unsigned i )
{
	// (ro, es, te) mode, here sum for sample is number of solved problems with time < te 
	// get best predict time for current point with different variants of time limits
	double point_cur_predict_time = 0.0, point_best_predict_time = HUGE_DOUBLE, 
		   cur_probability, point_best_time_limit = 0.0;
	unsigned cur_solved_in_time, point_best_solved_in_time = 0;
	vector<double> predict_times;
	predict_times.resize( predict_time_limites.size() );
	
	/*string temp_file_name = "temp_file";
	ofstream ofile( temp_file_name, std::ios::out );
	ofile << "time_limit predict solved_in_time probability cur_cnf_in_set_count" << std::endl;*/
	
	//unsigned index = 0, point_best_index = 0;
	for ( auto &cur_time_limit : predict_time_limites ) {
		//ofile << cur_time_limit << " ";
		cur_solved_in_time = 0;
		for ( unsigned j = set_index_arr[i]; j < set_index_arr[i + 1]; j++ )
			if ( ( cnf_issat_arr[j] ) && ( cnf_real_time_arr[j] > 0.0 ) && ( cnf_real_time_arr[j] <= cur_time_limit ) )
				cur_solved_in_time++;
		if ( !cur_solved_in_time )
			continue;
		cur_probability = (double)cur_solved_in_time / (double)cur_cnf_in_set_count;
		point_cur_predict_time = pow( 2.0, (double)cur_var_num ) * cur_time_limit * 3.0 / cur_probability;
		if ( point_cur_predict_time < point_best_predict_time ) {
			point_best_predict_time   = point_cur_predict_time;
			point_best_solved_in_time = cur_solved_in_time;
			point_best_time_limit     = cur_time_limit;
			//point_best_index          = index;
		}
		//index++;
		/*ofile << point_cur_predict_time << " " << cur_solved_in_time << " " << 
			     cur_probability << " " << cur_cnf_in_set_count << std::endl;*/
	}
	
	//ofile.close();
	
	// increase range if best value in right limit 
	/*double last_value = predict_time_limites[predict_time_limites.size()-1];
	if ( point_best_index == predict_time_limites.size() - 1 ) {
		predict_time_limites.push_back( last_value + predict_time_limit_step );
		te = predict_time_limites[predict_time_limites.size()-1];
		ofstream deep_predict_file( deep_predict_file_name.c_str(), ios_base::out | ios_base::app );
		deep_predict_file << std::endl << "predict_time_limites updated" <<std::endl;
		for ( auto &x : predict_time_limites )
			deep_predict_file << x << " ";
		deep_predict_file << "te " << te << std::endl;
		deep_predict_file.close();
	}*/
	
	solved_in_time_arr[i] = point_best_solved_in_time;
	time_limit_arr[i]     = point_best_time_limit;
	
	return point_best_predict_time;
}		

//---------------------------------------------------------
bool MPI_Predicter :: GetPredict()
{
// set_status_arr == 0, if there are some unsolved sat-problem in set 
// set_status_arr == 1, if at least one sat-problem in set was STOPPED
// set_status_arr == 2, if all CNF in set are solved and it's BKV
// set_status_arr == 4, if all CNF in set are solved, not BKV, but control process couldn't catch to stop it
	if ( verbosity > 2 )
		std::cout << "GetPredict()" << std::endl;
	
	int solved_tasks_count = 0;
	unsigned cur_var_num, solved_in_sample_count, 
			 cur_cnf_to_skip_count = 0, 
			 cur_cnf_to_stop_count = 0;
	double cur_predict_time  = 0.0,
		   cur_sum_part_time = 0.0,
		   cur_med_part_time = 0.0,
		   time1 = 0,
		   cur_cnf_time = 0.0;
	unsigned long long temp_llint = 0;
	int set_index_bound = 0,
		cur_cnf_in_set_count = 0;
	double current_time = Minisat::cpuTime();
	
	cnf_to_stop_arr.clear(); // every time get stop-list again
	stringstream sstream;
	bool isTeBkvUpdated;
	
	// fill arrays of summary and median times in set of CNF
	for ( unsigned i = 0; i < decomp_set_arr.size(); i++ ) {
		if ( set_status_arr[i] > 0 )
			continue; // skip UNSAT, SAT and STOPPED
		cur_cnf_to_stop_count = 0; // every time create array again
		cur_cnf_to_skip_count = 0;
		// count of CNF in set == ( set_index_arr[i + 1] - set_index_arr[i] )
		cur_cnf_in_set_count = set_index_arr[i + 1] - set_index_arr[i];
		solved_in_sample_count = 0;
		
		// if sat-problem is being solved (status 0 or 1) get current time
		for ( unsigned j = set_index_arr[i]; j < set_index_arr[i + 1]; j++ ) {
			if ( cnf_status_arr[j] == 1 ) // if problem was stopped
				set_status_arr[i] = 1; // mark status STOP to set
			else if ( cnf_status_arr[j] > 1 )
				solved_in_sample_count++;
		} // for ( j = set_index_arr[i]; j < set_index_arr[i + 1]; j++ )
		
		if ( deep_predict )
			cur_var_num = decomp_set_arr[i].var_choose_order.size();
		else // if !deep_predict
			cur_var_num = predict_to - i;
		
		if ( ( cur_var_num <= MAX_STEP_CNF_IN_SET ) && ( solved_in_sample_count > (unsigned)( 1 << cur_var_num ) ) ) {
			cerr << "solved_in_sample_count > ( 1 << cur_var_num )" << endl;
			cerr << solved_in_sample_count << " > " << ( 1 << cur_var_num ) << endl;
			return false;
		}
		
		if ( solved_in_sample_count == cur_cnf_in_set_count ) // if all CNF in set were solved
			set_status_arr[i] = 4; // then mark status SOLVED to set
		
		// further will be only sets whisch are being solved right now
		solved_cnf_count_arr[i] = solved_in_sample_count;

		sum_time_arr[i] = 0.0; // init value of sum - every time must start from 0.0
		//max_real_time_sample = 0;

		isTeBkvUpdated = false;
		
		// get current predict time
		if ( te == 0.0 ) {
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

			med_time_arr[i] = sum_time_arr[i] / (double)cur_cnf_in_set_count;
			cur_predict_time = med_time_arr[i] / (double)proc_count;
			cur_predict_time *= pow( 2.0, (double)cur_var_num );
		}
		else if ( te > 0.0 ) // here med_time_arr in (0,1)
			cur_predict_time = getCurPredictTime( cur_var_num, cur_cnf_in_set_count, i );
		
		predict_time_arr[i] = cur_predict_time;
		
		// if in set all sat-problems solved and unsat then set can has best predict
		if ( ( set_status_arr[i] == 4  ) && 
			 ( predict_time_arr[i] > 0 ) &&
			 ( ( best_predict_time == 0.0 ) ||
			   ( ( best_predict_time > 0.0 ) && ( ( predict_time_arr[i] < best_predict_time ) ) ) || 
			   ( isTeBkvUpdated ) ||
			   ( checkSimulatedGranted( predict_time_arr[i] ) )
			 )
		   )
		{
			// if all sat-problems were solved then it can be best predict
			set_status_arr[i]     = 2;
			IsRecordUpdated       = true;
			best_predict_time     = predict_time_arr[i];
			best_sum_time         = sum_time_arr[i];
			best_solved_in_time   = solved_in_time_arr[i];
			best_cnf_in_set_count = cur_cnf_in_set_count;
			best_time_limit       = time_limit_arr[i];
			
			//prev_best_sum = sum_time_arr[i];
			//prev_best_decomp_set_power = decomp_set_arr[i].var_choose_order.size();
			
			if ( deep_predict ) // Write info about new point in deep mode
				NewRecordPoint( i );
			//if ( IsRestartNeeded ) // don't stop immidiatly, new record can be found in calculated points
			//	return true;
		}
		else if ( ( best_predict_time > 0.0 ) && ( predict_time_arr[i] >= best_predict_time  ) && 
			      ( predict_time_arr[i] != HUGE_DOUBLE ) && (!isSolverSystemCalling) ) { // stop, predict >= best
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

	if ( verbosity > 2 )
		std::cout << "In GetPredict() after main loop" << std::endl;
	
	//cout << "before cycle of AddNewUncheckedArea()" << endl;
	if ( deep_predict >= 5 ) {
		// add to L2 once for every set
		for ( unsigned i = 0; i < set_status_arr.size(); ++i ) {
			if ( ( decomp_set_arr[i].IsAddedToL2 == false ) && ( set_status_arr[i] > 0 ) )  {
				decomp_set_arr[i].IsAddedToL2 = true;
				if ( set_status_arr[i] == 1 )
					global_stopped_points_count++;
				else 
					global_checked_points_count++;
				boost::dynamic_bitset<> bs = IntVecToBitsetPredict( decomp_set_arr[i].var_choose_order );
				AddNewUncheckedArea( bs, sstream );
			}
		}
		fstream deep_predict_file( deep_predict_file_name.c_str(), ios_base::out | ios_base::app );
		deep_predict_file << sstream.rdbuf();
		deep_predict_file.close();
	}
	
	whole_get_predict_time += Minisat::cpuTime() - current_time;

	if ( verbosity > 2 )
		cout << "In GetPredict() after AddNewUncheckedArea()" << endl;

	return true;
}

//---------------------------------------------------------
bool MPI_Predicter :: WritePredictToFile( int all_skip_count, double whole_time_sec )
{
// Write info about predict to file
	ofstream predict_file;
	predict_file.open( predict_file_name.c_str( ), ios :: out ); // create and open for writing
	if ( !( predict_file.is_open( ) ) )
	{ std::cout << "Error in opening of predict_file " << predict_file_name << std::endl; return false; }
	
	stringstream sstream;
	sstream << "Predict from "             << predict_from << std::endl
		    << "Predict to "               << predict_to << std::endl
			<< "Processor count "          << proc_count << std::endl
			<< "Count of CNF in set "      << cnf_in_set_count << std::endl
			<< "Count of CNF in all sets " << all_tasks_count << std::endl
			<< "Count of core-variables "  << core_len << std::endl
			<< "Solver name "              << solver_name << std::endl
			<< "Schema type "              << schema_type << std::endl
			<< "start activity "           << start_activity << std::endl;

	if ( deep_predict ) {
		sstream << "deep_predict "     << deep_predict << std::endl;
		sstream << "decomp_set_arr.size() " << decomp_set_arr.size();	
	}
	
	predict_file << sstream.rdbuf( );
	
	double med_cnf_time, min_cnf_time, max_cnf_time;
	double sample_variance; // sample variance for estimation of derivation
	unsigned max_time_cnf_value_mask = 0;
	bool isAllSolved;
	unsigned solved_count;
	for ( unsigned i = 0; i < decomp_set_arr.size(); i++ ) {
		med_cnf_time = 0; 
		min_cnf_time = 0; 
		max_cnf_time = 0;
		sample_variance = 0;
		sstream.str( "" );
		sstream.clear();
		sstream << endl << " ";
		solved_count = 0;
		
		if ( deep_predict )
			sstream << decomp_set_arr[i].var_choose_order.size();
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
		
		med_cnf_time = 0.0;
		// prepare start min and max values
		bool IsFirstNonNullFinded = false;
		unsigned count0 = 0, count1 = 0, count2 = 0, count3 = 0, count4 = 0, count5 = 0, count6 = 0, count7 = 0, count8 = 0;
		isAllSolved = true;
		min_cnf_time = max_cnf_time = -1.0;
		for ( unsigned j = set_index_arr[i]; j < set_index_arr[i + 1]; ++j ) {
			if ( cnf_status_arr[j] <= 1 ) { // skip unsolved and stopped problems
				isAllSolved = false;
				continue;
			}
			else
				solved_count++;
			if ( ( !IsFirstNonNullFinded ) && ( cnf_real_time_arr[j] > 0.0 ) ) {
				min_cnf_time = max_cnf_time = cnf_real_time_arr[j];
				IsFirstNonNullFinded = true;
			}
			med_cnf_time += cnf_real_time_arr[j];
			if ( cnf_real_time_arr[j] < min_cnf_time )
				min_cnf_time = cnf_real_time_arr[j];
			if ( cnf_real_time_arr[j] > max_cnf_time )
				max_cnf_time = cnf_real_time_arr[j];
			if ( cnf_prepr_arr[j] )                  count0++;
			else if ( cnf_real_time_arr[j] < 0.001 ) count1++;
			else if ( cnf_real_time_arr[j] < 0.01 )  count2++;
			else if ( cnf_real_time_arr[j] < 0.1  )  count3++;
			else if ( cnf_real_time_arr[j] < 1    )  count4++;
			else if ( cnf_real_time_arr[j] < 10   )  count5++;
			else if ( cnf_real_time_arr[j] < 100  )  count6++;
			else if ( cnf_real_time_arr[j] < 1000 )  count7++;
			else count8++;
		}
		
		/*if ( ( set_status_arr[i] > 0 ) && ( solved_cnf_count_arr[i] != solved_count ) ) {
			cerr << "solved_cnf_count_arr[i] != solved_count" << endl;
			cerr << solved_cnf_count_arr[i] << " != " << solved_count << endl;
			cerr << "set_status_arr[i] " << set_status_arr[i] << endl;
			cerr << "cnf_status_arr" << endl;
			for ( unsigned j = set_index_arr[i]; j < set_index_arr[i + 1]; ++j )
				cerr << cnf_status_arr[j] << endl;
			exit(1);
		}*/
		if ( (isAllSolved) && (( solved_cnf_count_arr[i] > 0 )) )
			med_cnf_time /= solved_cnf_count_arr[i];
		else
			med_cnf_time = -1.0;
		
		// compute sample_variancse
		if ( ( med_cnf_time > 0.0 ) && ( solved_cnf_count_arr[i] > 0 ) ) {
			for ( unsigned j = set_index_arr[i]; j < set_index_arr[i + 1]; j++ )
				sample_variance += pow(cnf_real_time_arr[j] - med_cnf_time, 2);
			sample_variance /= ( solved_cnf_count_arr[i] - 1);
		}
		else
			sample_variance = -1.0;
		
		if ( sample_variance > start_sample_variance_limit ) {
			/*cout << "new sample_varianse" << endl;
			cout << "sample_variance " << sample_variance << endl;
			cout << "var_choose_order.size() " << var_choose_order.size() << endl;*/
			start_sample_variance_limit *= 2;
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
		
		sstream << " prepr: "       << count0;
		sstream << " (0, 0.001): "  << count1;
		sstream << " (0.001, 0.01): " << count2;
		sstream << " (0.01, 0.1): " << count3;
		sstream << " (0.1, 1): "    << count4;
		sstream << " (1, 10): "     << count5;
		sstream << " (10, 100): "   << count6;
		sstream << " (100, 1000): " << count7;
		sstream << " (1000, .): "   << count8;

		predict_file << sstream.rdbuf( );
	}

	sstream.str( "" ); // clear stringstream
	sstream.clear();
	sstream <<std::endl << std::endl;
	sstream << "All skipped count "       << all_skip_count << std::endl;
	sstream << "Current solved tasks count " << solved_tasks_count << std::endl;
	sstream << "Best var num "               << best_var_num << std::endl;
	sstream << "Best predict time " << left << scientific << best_predict_time << std::endl;
	sstream << "Predict took time "          << whole_time_sec << std::endl;
	predict_file << sstream.rdbuf( );
	predict_file.close( );
	
	return true;
}

/*
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
*/

void MPI_Predicter :: AllocatePredictArrays()
{
	all_tasks_count = 0;
	unsigned uint;
	int cur_tasks_count;
	// array of random set lengths 
	set_len_arr.clear();
	for ( unsigned i = 0; i < decomp_set_arr.size(); i++ ) {
		if ( decomp_set_arr[i].var_choose_order.size() > MAX_STEP_CNF_IN_SET ) // if too many vars
			cur_tasks_count = cnf_in_set_count;
		else {
			uint = ( 1 << decomp_set_arr[i].var_choose_order.size() ); // current count of all tasks
			cur_tasks_count = ( cnf_in_set_count < (int)uint ) ? cnf_in_set_count : uint;
		}
	    set_len_arr.push_back( cur_tasks_count );
		all_tasks_count += cur_tasks_count;
	}
	/*ofstream ofile("set_len_arr", ios_base::out);
	for ( unsigned i=0; i < set_len_arr.size(); ++i )
		ofile << i << " " << set_len_arr[i] << endl;
	ofile.close();*/
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
	double current_time = Minisat::cpuTime();

	if ( verbosity > 1 ) {
		cout << "Started AddNewUncheckedArea() with" << endl;
		//sstream << point.to_string() << endl;
	}

	new_ua.center = point;
	new_ua.radius = 1;
	new_ua.checked_points.resize( point.size() );
	
	for ( L1_it = L1.begin(); L1_it != L1.end(); L1_it++ ) {
		// necessary condition - points with close count of 1s
		if ( abs( (int)(*L1_it).center.count() - (int)new_ua.center.count() ) <= (*L1_it).radius ) { 
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
		if ( abs( (int)(*L2_it).center.count() - (int)new_ua.center.count() ) <= (*L2_it).radius ) { 
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
	
	//to_string( new_ua.checked_points, str );
	//sstream << str << endl;
	if ( ( new_ua.checked_points.count() == 0 ) && ( !isFirstPoint ) )
		sstream << "***Error. new_ua.center.count() == 0" << endl;
	if ( verbosity > 1 )
		cout << "In AddNewUncheckedArea() end" << endl;
	
	whole_add_new_unchecked_area_time += Minisat::cpuTime() - current_time;
}

bool MPI_Predicter :: GetDeepPredictTasks( )
{
// Make tasks for checking neighbours of current point
	// var_choose_order - current best decomosition (point)
	cout << endl << "*** GetDeepPredictTasks" << endl;
	cout << "global_deep_point_index " << global_deep_point_index << endl;

	unsigned points_to_check;

	if ( isFirstPoint )// compute value of function only in first point
		points_to_check = total_decomp_set_count = 1;
	else {
		// core_len variants of vectors with Hamming distanse 
		// update - make points only for current vars count 
		if ( global_deep_point_index == 0 ) { // calculate only at 1st time
			// get combinations to constuct all vectors on given Hamming distanse
			MakeCombinations( core_len, cur_vars_changing, combinations );
			random_shuffle( combinations.begin(), combinations.end() ); // shuffle vars for getting equal random of changing
			total_decomp_set_count = combinations.size();
			cout << "total_decomp_set_count " << total_decomp_set_count << endl;
		}
		points_to_check = total_decomp_set_count - global_deep_point_index; // how many points to check
	}
	
	prev_area_best_predict_time = best_predict_time;
	unsigned cur_index;
	boost::dynamic_bitset<> new_point;

	//if ( ( deep_predict == 6 ) && ( !isFirstPoint ) )
	//	deep_predict_cur_var = current_unchecked_area.center.count();
	
	if ( verbosity > 1 ) {
		cout << "AllocateDeepArrays() done" << endl;
		cout << "decomp_set_arr.size() " << decomp_set_arr.size() << endl;
	}
	decomp_set_arr.clear();
	decomp_set d_s;
	unsigned current_skipped = 0;
	stringstream sstream;
	vector<int> new_var_choose_order;
	//sstream << "current point to check" << endl;
	for ( unsigned i = 0; i < points_to_check; ++i ) { // several decomp sets
		if ( isFirstPoint )
			new_var_choose_order = var_choose_order;
		else {
			cur_index = global_deep_point_index + i;
			/*if ( deep_predict <= 2 )
				ChangeVarChooseOrder( var_choose_order, cur_vars_changing, cur_var_count, new_var_choose_order );
			else if ( deep_predict != 6 )
				GetNewHammingPoint( var_choose_order, cur_vars_changing, cur_var_count, combinations[cur_index],
									new_var_choose_order );*/
			if ( deep_predict >= 5 ) {
				if ( current_unchecked_area.checked_points[cur_index] == 1 ) {
					current_skipped++;
					continue; // checked already
				}
				
				new_point = current_unchecked_area.center; // new point based on center point
				new_point[cur_index] = ( current_unchecked_area.center[cur_index] == 1 ) ? 0 : 1;
				new_var_choose_order = BitsetToIntVecPredict( new_point );
			}
		}
		
		d_s.var_choose_order = new_var_choose_order;
		d_s.cur_var_changing = cur_vars_changing;
		d_s.IsAddedToL2      = false;
		decomp_set_arr.push_back( d_s ); // add new decomp set
	} // for ( int i = 0; i < decomp_set_count; i++ )
	
	if ( ts_strategy == 0 )
		random_shuffle( decomp_set_arr.begin(), decomp_set_arr.end() );
	else if ( ts_strategy == 1 ) {
		// sort decomp sets by activity
		vector<decomp_set> :: iterator dec_it;
		vector<int> :: iterator vec_it;
		for ( dec_it = decomp_set_arr.begin(); dec_it != decomp_set_arr.end(); ++dec_it ) {
			(*dec_it).med_var_activity = 0;
			for ( vec_it = (*dec_it).var_choose_order.begin(); vec_it != (*dec_it).var_choose_order.end(); ++vec_it )
				(*dec_it).med_var_activity += total_var_activity[ core_var_indexes.find(*vec_it)->second ];
			(*dec_it).med_var_activity /= (*dec_it).var_choose_order.size();
		}
		sort( decomp_set_arr.begin(), decomp_set_arr.end(), ds_compareByActivity );
		if ( verbosity > 0 ) {
			cout << "decomp_set_arr activity" << endl;
			for ( dec_it = decomp_set_arr.begin(); dec_it != decomp_set_arr.end(); ++dec_it )
				cout << (*dec_it).med_var_activity << endl;
		}
	}

	global_skipped_points_count += current_skipped;

	// Common for all procedures of getting deep tasks
	AllocatePredictArrays();

	if ( verbosity > 0 )
		cout << "AllocatePredictArrays() done " << endl;
	
	sstream << "total_decomp_set_count " << total_decomp_set_count << endl;
	sstream << "global_deep_point_index " << global_deep_point_index << endl;
	sstream << "must be checked " << decomp_set_arr.size() << " from " << points_to_check << endl;
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
	sstream << "whole_deep_time "             << Minisat::cpuTime() - whole_deep_time << " s" << endl;
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
	
	if ( verbosity > 0 )
		cout << sstream.str();
	
	global_deep_point_index += decomp_set_arr.size() + current_skipped;
	
	sstream << endl;
	fstream deep_predict_file( deep_predict_file_name.c_str(), ios_base::out | ios_base::app );
	deep_predict_file << sstream.rdbuf();
	deep_predict_file.close();

	return true;
}