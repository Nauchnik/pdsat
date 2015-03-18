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
	max_solving_time_koef       ( 0 ),
	finding_first_sat_time      ( 0 ),
	total_start_time            ( 0 ),
	no_increm                   ( false )
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
// Solve subproblems with MPI
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &corecount );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	isPredict = false;
	double iteration_start_time, iteration_final_time;
	double whole_final_time;
	std::stringstream sstream;
	int break_message = -1;
	int stop_message  = -2;
	
	if ( corecount < 2 ) { 
		std::cerr << "corecount < 2"; 
		MPI_Abort( MPI_COMM_WORLD, 0 );
	}
	
	if ( rank != 0 ) { // computing processes
		if ( !ComputeProcessSolve() ) {
			std::cerr << "in ComputeProcessSovle" << std::endl; 
			MPI_Abort( MPI_COMM_WORLD, 0 );
		}
	}
	else { // rank == 0, control process
		std::cout << "*** MPI_Solve is running ***" << std::endl;
		total_start_time = MPI_Wtime();
		//for (;;) {
			sstream << base_solving_info_file_name << "_" << solving_iteration_count;
			solving_info_file_name = sstream.str();
			sstream.clear(); sstream.str("");
			std::cout << "solving_info_file_name " << solving_info_file_name << std::endl; 
			
			iteration_start_time = MPI_Wtime();
			
			std::vector<std::vector<bool>> interrupted_problems_var_values;
			std::vector< satisfying_assignment > satisfying_assignments;
			ControlProcessSolve( var_choose_order, interrupted_problems_var_values, satisfying_assignments );
			std::cout << "final interrupted_problems_var_values.size() " << interrupted_problems_var_values.size() << std::endl;
			std::cout << "final satisfying_assignments.size() " << satisfying_assignments.size() << std::endl;
			
			if ( interrupted_problems_var_values.size() > 0 ) {
				std::ofstream ofile("interrupted_problems");
				for ( auto &x : interrupted_problems_var_values ) {
					for ( unsigned j=0; j < x.size(); j++ )
						ofile << x[j];
					ofile << std::endl;
				}
				ofile.close();
			}

			unsigned ones_count;
			std::string tmp_str;
			if ( satisfying_assignments.size() > 0 ) {
				std::ofstream ofile("satisfying_assignments");
				for ( auto &x : satisfying_assignments )
					ofile << x.str << std::endl;
				ofile.close(); ofile.clear();
				ofile.open( "decomp_set_satisfying_assignments" );
				for ( auto &x : satisfying_assignments ) {
					ones_count = 0;
					tmp_str = "";
					for ( auto &y : var_choose_order ) {
						if ( x.str[y-1] == '1' ) 
							ones_count++;
						tmp_str += x.str[y-1];
					}
					ofile << x.solving_time << " s " << "ones " << ones_count << " " << tmp_str << std::endl;
				}
				ofile.close();
			}
			
			iteration_final_time = MPI_Wtime() - iteration_start_time;
			WriteTimeToFile( iteration_final_time );
			
			/*if ( max_solving_time_koef > 0.0 ) { // if increasing time limit for same subproblems
				if ( ( sat_count && !IsSolveAll ) || ( max_solving_time_koef == 0.0 ) )
					break; // exit if SAT set found or iterative solving is not needed
				
				// send messages for breaking low loop on compute processes
				if ( interrupted_count ) {
					std::cout << "sending break messages to all computing processes" << std::endl;
					for ( int i = 1; i < corecount; i++ )
						MPI_Send( &break_message, 1, MPI_INT, i, 0, MPI_COMM_WORLD );
				}
				else 
					break;
				
				max_solving_time *= max_solving_time_koef; // increase time limit
			}
			else break;*/
		//}
		
		whole_final_time = MPI_Wtime() - total_start_time;
		solving_info_file_name = base_solving_info_file_name + "_total";
		WriteTimeToFile( whole_final_time );
		
		// send messages for finalizing
		for ( int i = 1; i < corecount; i++ )
			MPI_Send( &stop_message, 1, MPI_INT, i, 0, MPI_COMM_WORLD );
		
		MPI_Finalize( );
	}

	return 0;
}

void MPI_Solver :: AddSolvingTimeToArray( ProblemStates cur_problem_state, double cnf_time_from_node, 
										  double *solving_times )
{
	// solving_times[0]  == min
	// solving_times[1]  == max
	// solving_times[2]  == med
	// solving_times[3]  == sat
	if( cur_problem_state != Interrupted ){
		if ( cnf_time_from_node < 0 ) {
			std::cerr << "cnf_time_from_node " << cnf_time_from_node << " < 0" << std::endl;
			exit(1);
		}
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

bool MPI_Solver :: SolverRun( Solver *&S, unsigned long long &process_sat_count, 
							  double &cnf_time_from_node, int current_task_index,
							  std::vector< std::vector<bool> > &interrupted_problems_var_values_from_process,
							  std::vector< std::vector<bool> > &sat_assignment_from_process )
{
	if ( verbosity > 1 )
		std::cout << "start SolverRun()" << std::endl;
	
	interrupted_problems_var_values_from_process.clear();
	sat_assignment_from_process.clear();
	process_sat_count = 0;
	vec< vec<Lit> > dummy_vec;
	lbool ret;
	double total_time = 0;
	unsigned long long current_tasks_solved = 0;
	ProblemStates cur_problem_state;
	std::stringstream sstream;
	std::ofstream ofile;
	unsigned long long batch_interrupted_count = 0;

	solving_times[0] = 1 << 30; // start min len
	for ( unsigned i = 1; i < SOLVING_TIME_LEN; i++ )
		solving_times[i] = 0;

	//vector< boost::dynamic_bitset<> > vec_bitset;
	//boost::dynamic_bitset<> cur_bitset;
	std::vector<bool> cur_interrupted_values;
	cur_interrupted_values.resize( var_choose_order.size() );
	
	unsigned long long before_binary_length = 0;
	
	MakeAssignsFromMasks( full_mask, part_mask, mask_value, dummy_vec );
	
	if ( ( verbosity > 1 ) && ( rank == 1 ) ) {
		std::cout << "dummy_vec size" << dummy_vec.size() << std::endl;
		for ( int i = 0; i < dummy_vec.size(); ++i ) {
			for ( int j=0; j < dummy_vec[i].size(); ++j )
				std::cout << dummy_vec[i][j].x << " ";
			std::cout << std::endl;
		}
	}
	
	for ( int i = 0; i < dummy_vec.size(); ++i )
		if ( dummy_vec[i].size() == 0 ) {
			std::cerr << "dummy_vec.size() == 0" << std::endl;
			return false;
		}
	
	uint64_t prev_starts, prev_conflicts;
	for ( int i=0; i < dummy_vec.size(); ++i ) {
		// save current state to check differences
		prev_starts    = S->starts;
		prev_conflicts = S->conflicts;
			
		cnf_time_from_node = Minisat :: cpuTime();
		ret = S->solveLimited( dummy_vec[i] );
		cnf_time_from_node = Minisat :: cpuTime() - cnf_time_from_node;

		if ( no_increm )
			S->clearDB(); // clear database if incremental solving disabled
			
		//ret = S->solveLimited( dummy_vec[i], true, false ); // for SimpSolver

		total_time += cnf_time_from_node;
			
		if ( ret == l_Undef )
			cur_problem_state = Interrupted; // interrupted cause of restarts or time limit
		else if ( ( S->starts - prev_starts <= 1 ) && ( S->conflicts == prev_conflicts ) )
			cur_problem_state = SolvedOnPreprocessing;  // solved by BCP
		else
			cur_problem_state = Solved; // just solved
		
		AddSolvingTimeToArray( cur_problem_state, cnf_time_from_node, solving_times );
		
		if ( cur_problem_state == Interrupted ) {
			batch_interrupted_count++;
			if ( cur_interrupted_values.size() < (unsigned)dummy_vec[i].size() ) {
					std::cerr << "cur_interrupted_values.size() < dummy_vec[i].size()" << std::endl;
					std::cerr << cur_interrupted_values.size() << " < " << dummy_vec[i].size() << std::endl;
					return false;
				}
				for ( int j = 0; j < dummy_vec[i].size(); ++j )
					cur_interrupted_values[j] = ( dummy_vec[i][j].x % 2 == 0 ) ? true : false; 
				interrupted_problems_var_values_from_process.push_back( cur_interrupted_values );
			}
			
		if ( ret == l_True ) {
			process_sat_count++;
			std::cout << "process_sat_count " << process_sat_count << std::endl;
			if ( !solving_times[3] ) {
				solving_times[3] = cnf_time_from_node; // time of 1st SAT, write only once
				std::cout << "SAT time " << solving_times[3] << std::endl;
			}
			b_SAT_set_array.resize( S->model.size() );
			for ( int i=0; i < S->model.size(); i++ )
				b_SAT_set_array[i] = ( S->model[i] == l_True) ? true : false;
			// check res file for SAT set existing
			if ( !AnalyzeSATset( cnf_time_from_node ) ) {
				// is't needed to deallocate memory - MPI_Abort will do it	
				std::cout << "Error in Analyzer" << std::endl;
				MPI_Abort( MPI_COMM_WORLD, 0 );
				return false;
			}
			if ( !IsSolveAll )
				break;
			sat_assignment_from_process.push_back( b_SAT_set_array );
			std::cout << "satisfying assignment.size() " << b_SAT_set_array.size() << std::endl;
			for ( unsigned j=0; j < b_SAT_set_array.size() ; j++ )
				std::cout << b_SAT_set_array[j];
		}
	}

	/*if ( batch_interrupted_count ) {
		std::fstream file( new_assumptions_file_name.c_str(), std::ios_base::in );
		if ( file.peek() == std::fstream::traits_type::eof() ) { // if file is empty
			file.close();
			file.open( new_assumptions_file_name.c_str(), std::ios_base::out );
			file << var_choose_order.size();  // write length of boolean vectors
			file.close();
		}
		else file.close();
			
		ofile.open( new_assumptions_file_name.c_str(), std::ios_base::out | std::ios_base::app | std::ios_base::binary );
		unsigned long long ul;
		for( unsigned i=0; i < vec_bitset.size(); ++i ) {
			//ul = vec_bitset[i].to_ullong();
			ul = BitsetToUllong( vec_bitset[i] );
			ofile.write( (char*)&ul, sizeof(ul) );
		}
		ofile.close();
	}*/
	
	solving_times[2] = total_time / dummy_vec.size(); // med time for current batch

	return true;
}

void MPI_Solver :: WriteSolvingTimeInfo( double *solving_times, unsigned solved_tasks_count )
{
	int time_sec, time_minutes, time_hours;
	
	if ( solving_times[0] < 0 ) {
		std::cerr << "solving_times[0] < 0" << std::endl;
		exit(1);
	}
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

	std::stringstream sstream;
	sstream << "no_increm " << no_increm << std::endl;
	sstream << "var_choose_order size " << var_choose_order.size() << std::endl;
	for ( unsigned i=0; i < var_choose_order.size(); i++ )
		sstream << var_choose_order[i] << " ";

	sstream << std::endl;
	sstream << "solved_tasks_count     " << solved_tasks_count    << " / " << all_tasks_count << std::endl;
	sstream << "solved_problems_count  " << solved_problems_count << std::endl;
	sstream << "*** total_solving_times" << std::endl;
	sstream << "min                    " << total_solving_times[0] << std::endl;
	sstream << "max                    " << total_solving_times[1] << std::endl;
	sstream << "med                    " << total_solving_times[2] << std::endl;
	sstream << "sat                    " << total_solving_times[3] << std::endl;
	sstream << "finding_first_sat_time ";
	cpuTimeInHours( finding_first_sat_time, time_hours, time_minutes, time_sec );
	sstream << time_hours << " h " << time_minutes << " m " << time_sec << " s" << std::endl;
	sstream << "sat_count              " << sat_count << std::endl;
	if ( solved_problems_count > 0 ) {
		sstream << "SOLVED BY PREPROCESSING, count ";
		sstream.width( 12 ); 
		sstream << std::left << std::fixed << (unsigned)total_solving_times[4]  << " " << total_solving_times[4] * 100 / solved_problems_count << " %" << std::endl;
		sstream << "(0,   0.0001)   count ";
		sstream.width( 12 ); 
		sstream << std::left << std::fixed << (unsigned)total_solving_times[5]  << " " << total_solving_times[5] * 100 / solved_problems_count << " %" << std::endl;
		sstream << "(0.0001, 0.001) count ";
		sstream.width( 12 ); 
		sstream << std::left << std::fixed << (unsigned)total_solving_times[6]  << " " << total_solving_times[6] * 100 / solved_problems_count << " %" << std::endl;
		sstream << "(0.001, 0.01)   count ";
		sstream.width( 12 );
		sstream << std::left << std::fixed << (unsigned)total_solving_times[7]  << " " << total_solving_times[7] * 100 / solved_problems_count << " %" << std::endl;		
		sstream << "(0.01, 0.1)     count ";
		sstream.width( 12 ); 
		sstream << std::left << std::fixed << (unsigned)total_solving_times[8]  << " " << total_solving_times[8] * 100 / solved_problems_count << " %" << std::endl;
		sstream << "(0.1,    1)     count ";
		sstream.width( 12 ); 
		sstream << std::left << std::fixed << (unsigned)total_solving_times[9]  << " " << total_solving_times[9] * 100 / solved_problems_count << " %" << std::endl;
		sstream << "(1,     10)     count "; 
		sstream.width( 12 ); 
		sstream << std::left << std::fixed << (unsigned)total_solving_times[10] << " " << total_solving_times[10] * 100 / solved_problems_count << " %" << std::endl;
		sstream << "(10,   100)     count "; 
		sstream.width( 12 ); 
		sstream << std::left << std::fixed << (unsigned)total_solving_times[11] << " " << total_solving_times[11] * 100 / solved_problems_count << " %" << std::endl;
		sstream << "(100, 1000)     count ";
		sstream.width( 12 );
		sstream << std::left << std::fixed << (unsigned)total_solving_times[12] << " " << total_solving_times[12] * 100 / solved_problems_count << " %" << std::endl;
		sstream << "(1000, inf)     count "; 
		sstream.width( 12 );
		sstream << std::left << std::fixed << (unsigned)total_solving_times[13] << " " << total_solving_times[13] * 100 / solved_problems_count << " %" << std::endl;
		sstream << "INTERRUPTED, solving time > " << max_solving_time << " ";
		sstream.width( 12 );
		sstream << std::left << std::fixed << (unsigned)total_solving_times[14] << " " << total_solving_times[14] * 100 / solved_problems_count << " %" << std::endl;		
	}
	if ( verbosity > 0 )
		std::cout << sstream.str();
	std::ofstream solving_info_file( solving_info_file_name.c_str() );
	solving_info_file << sstream.str();
	solving_info_file.close();
	sstream.clear(); sstream.str("");
}

bool MPI_Solver :: ControlProcessSolve( std::vector<int> extern_var_choose_order, 
									    std::vector<std::vector<bool>> &interrupted_problems_var_values,
										std::vector<satisfying_assignment> &satisfying_assignments )
{
	interrupted_problems_var_values.clear();
	std::vector<bool> cur_interrupted_problems_var_values;
	satisfying_assignment cur_satisfying_assignment;
	std::cout << std::endl << "ControlProcessSolve is running" << std::endl;
	std::cout << "solving_iteration_count " << solving_iteration_count << std::endl;
	
	if ( solving_iteration_count == 0 ) {
		if ( !ReadIntCNF() ) { // Read original CNF
			std::cerr << "Error in ReadIntCNF" << std::endl; 
			MPI_Abort( MPI_COMM_WORLD, 0 );
		}
		if ( !MakeVarChoose() ) { 
			std::cerr << "Error in MakeVarChoose" << std::endl; 
			MPI_Abort( MPI_COMM_WORLD, 0 );
		}
	}

	if ( extern_var_choose_order.size() > 0 )
		var_choose_order = extern_var_choose_order;

	std::cout << "var_choose_order " << std::endl;
	for ( auto &x : var_choose_order )
		std::cout << x << " ";
	std::cout << std::endl;
	
	// log(a)/log(b) = log(_a)b
	unsigned max_possible_tasks_count = (unsigned)(pow( 2, ceil( log(corecount - 1)/log(2) ))) * (unsigned)(pow(2,koef_val) );
	std::cout << "max_possible_tasks_count " << max_possible_tasks_count << std::endl;
	std::cout << "current part_mask_var_count " << part_mask_var_count << std::endl; 
	part_mask_var_count = (unsigned)(log(max_possible_tasks_count)/log(2));
	if ( part_mask_var_count > var_choose_order.size() )
		part_mask_var_count = var_choose_order.size();
	
	// change batch size to treshold value if needed
	if ( var_choose_order.size() - part_mask_var_count > RECOMMEND_BATCH_VAR_COUNT ) {
		part_mask_var_count = var_choose_order.size() - RECOMMEND_BATCH_VAR_COUNT;
		std::cout << "part_mask_var_count changed to " << part_mask_var_count << std::endl;
	}
	if ( part_mask_var_count > MAX_PART_MASK_VAR_COUNT )
		part_mask_var_count = MAX_PART_MASK_VAR_COUNT;
	if ( var_choose_order.size() - part_mask_var_count > MAX_BATCH_VAR_COUNT ) {
		std::cerr << "Error. var_choose_order.size() - part_mask_var_count > MAX_BATCH_VAR_COUNT" << std::endl;
		std::cerr << var_choose_order.size() - part_mask_var_count << " < " << MAX_BATCH_VAR_COUNT << std::endl;
		MPI_Abort( MPI_COMM_WORLD, 0 );
	}
	std::cout << "part_mask_var_count " << part_mask_var_count << std::endl;

	// get default count of tasks = power of part_mask_var_count
	unsigned part_var_power = ( 1 << part_mask_var_count );
	std::cout << "part_var_power " << part_var_power << std::endl;
	// TODO add extended tasks counting
	all_tasks_count = part_var_power;
	std::cout << "all_tasks_count " << all_tasks_count << std::endl;
	
	if ( (int)all_tasks_count < corecount-1 ) {
		std::cerr << "Error. all_tasks_count < corecount-1" << std::endl; 
		std::cerr << all_tasks_count << " < " << corecount-1 << std::endl;
		MPI_Abort( MPI_COMM_WORLD, 0 );
	}
	
	if ( solving_iteration_count == 0 )
		PrintParams( );
	
	if ( skip_tasks >= all_tasks_count ) {
		std::cerr << "skip_tasks >= all_tasks_count " << std::endl;
		std::cerr << skip_tasks << " >= " << all_tasks_count << std::endl;
		MPI_Abort( MPI_COMM_WORLD, 0 );
	}
	
	values_arr.resize( all_tasks_count );
	for ( unsigned i = 0; i < values_arr.size(); ++i )
		values_arr[i].resize( FULL_MASK_LEN );
	
	if ( !MakeStandardMasks( part_var_power ) ) {
		std::cerr << "Error in MakeStandartMasks" << std::endl; 
		MPI_Abort( MPI_COMM_WORLD, 0 );
	}
	std::cout << "Correct end of MakeStandartMasks" << std::endl;
	
	unsigned solved_tasks_count = 0;
	// write init info
	WriteSolvingTimeInfo( solving_times, solved_tasks_count );
	int *var_choose_order_int = new int[MAX_CORE_LEN];
	for( unsigned i=0; i < MAX_CORE_LEN; ++i ) {
		if ( i < var_choose_order.size() )
			var_choose_order_int[i] = var_choose_order[i];
		else 
			var_choose_order_int[i] = -1;
	}
	
	std::cout << "before sending configuration info" << std::endl;
	// send core_len once to every compute process
	for ( int i=0; i < corecount-1; ++i ) {
		MPI_Send( &core_len,                1, MPI_INT,                 i + 1, 0, MPI_COMM_WORLD );
		MPI_Send( &all_tasks_count,         1, MPI_UNSIGNED,            i + 1, 0, MPI_COMM_WORLD );
		MPI_Send( &solving_iteration_count, 1, MPI_INT,                 i + 1, 0, MPI_COMM_WORLD );
		MPI_Send( &max_solving_time,        1, MPI_DOUBLE,              i + 1, 0, MPI_COMM_WORLD );
		MPI_Send( &start_activity,          1, MPI_DOUBLE,              i + 1, 0, MPI_COMM_WORLD );
		MPI_Send( var_choose_order_int,     MAX_CORE_LEN, MPI_INT, i + 1, 0, MPI_COMM_WORLD );
	}
	delete[] var_choose_order_int;
	
	int next_task_index = 0;
	
	//unsigned start_tasks_count = min( corecount - 1, (int)all_tasks_count );
	unsigned start_tasks_count = ((corecount - 1) < (int)all_tasks_count) ? (unsigned)(corecount - 1) : all_tasks_count;
	int char_arr_len;
	char *char_arr;
	unsigned elem_index;
	
	std::cout << "start_tasks_count " << start_tasks_count << std::endl;
	// send to all cores (except # 0) tasks from 1st range
	for ( int i = 0; i < (int)start_tasks_count; ++i ) {
		// send new index of task for reading tasks from file
		MPI_Send( &next_task_index, 1, MPI_INT, i + 1, 0, MPI_COMM_WORLD );

		if ( ( verbosity > 1 ) && ( i == 0 ) )
			std::cout << "sended next_task_index " << next_task_index << std::endl;
		
		copy( values_arr[i].begin(), values_arr[i].end(), mask_value );
		MPI_Send( full_mask,  FULL_MASK_LEN, MPI_UNSIGNED, i + 1, 0, MPI_COMM_WORLD );
		MPI_Send( part_mask,  FULL_MASK_LEN, MPI_UNSIGNED, i + 1, 0, MPI_COMM_WORLD );
		MPI_Send( mask_value, FULL_MASK_LEN, MPI_UNSIGNED, i + 1, 0, MPI_COMM_WORLD );
		next_task_index++;
	}
	std::cout << "after sending start_tasks_count" << std::endl;
	
	total_solving_times[0] = 1 << 30; // start min len
	for ( unsigned i = 1; i < total_solving_times.size(); ++i )
		total_solving_times[i] = 0;
	process_sat_count = 0;
	MPI_Status status, current_status;
	
	while ( solved_tasks_count < all_tasks_count ) {
		// recieve from core message about solved task 		
		MPI_Recv( &process_sat_count, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		if ( verbosity > 0 )
			std::cout << "recieved process_sat_count " << process_sat_count << std::endl;
		current_status = status;
		MPI_Recv( solving_times, SOLVING_TIME_LEN, MPI_DOUBLE, current_status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		if ( verbosity > 0 )
			std::cout << "recieved solving_times " << std::endl;
		
		// get interrupted tasks if such exist
		MPI_Probe( current_status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );;
		MPI_Get_count( &status, MPI_CHAR, &char_arr_len );
		if ( ( char_arr_len > 1 ) && ( char_arr_len % var_choose_order.size() != 0 ) ) {
			std::cerr << "char_arr_len % var_choose_order.size() != 0" << std::endl;
			MPI_Abort( MPI_COMM_WORLD, 0 );
		}
		if ( char_arr_len > 0 ) {
			char_arr = new char[char_arr_len];
			MPI_Recv( char_arr, char_arr_len, MPI_CHAR, current_status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
			if ( char_arr_len > 1 ) {
				//std::cout << "recieved char_arr_len " << char_arr_len << std::endl;
				cur_interrupted_problems_var_values.resize( var_choose_order.size() );
				elem_index=0;
				for ( int j=0; j < var_choose_order.size(); j++ ) {
					cur_interrupted_problems_var_values[elem_index++] = (char_arr[j] == '1' ? true : false);
					if ( (j+1) % var_choose_order.size() == 0 ) {
						interrupted_problems_var_values.push_back( cur_interrupted_problems_var_values );
						elem_index=0;
					}
				}
			}
			delete[] char_arr;
		}
		
		// get satisfying assignments if such exist
		MPI_Probe( current_status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );;
		MPI_Get_count( &status, MPI_CHAR, &char_arr_len );
		if ( ( char_arr_len > 1 ) && ( char_arr_len % var_count != 0 ) ) {
			std::cerr << "char_arr_len % var_count != 0" << std::endl;
			std::cerr << char_arr_len << " % " << var_count << " != 0 " << std::endl;
			MPI_Abort( MPI_COMM_WORLD, 0 );
		}
		if ( char_arr_len > 0 ) {
			cur_satisfying_assignment.solving_time = solving_times[3]; // sat solving time
			char_arr = new char[char_arr_len];
			MPI_Recv( char_arr, char_arr_len, MPI_CHAR, current_status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
			if ( char_arr_len > 1 ) { // read several assignments from one array
				std::cout << "recieved char_arr_len " << char_arr_len << std::endl;
				elem_index=0;
				cur_satisfying_assignment.str.resize( var_count );
				for ( int j=0; j < char_arr_len; j++ ) {
					cur_satisfying_assignment.str[elem_index++] = char_arr[j];
					if ( (j+1) % var_count == 0 ) {
						satisfying_assignments.push_back( cur_satisfying_assignment );
						elem_index=0;
					}
				}
			}
			delete[] char_arr;
		}
		
		if ( char_arr_len > 1 )
			std::cout << "interrupted_problems_var_values.size() " << interrupted_problems_var_values.size() << std::endl;

		/*char_send_array = new char[interrupted_problems_var_values_from_process.size() * var_choose_order.size()];
		char_send_array_index = 0;
		for ( auto &x : interrupted_problems_var_values_from_process )
			for ( auto &y : x )
				char_send_array[char_send_array_index++] = (y == true ? '1' : '0');
		MPI_Send( char_send_array, char_send_array_index, MPI_CHAR, 0, 0, MPI_COMM_WORLD );
		delete[] char_send_array;*/
		solved_tasks_count++;
		
		if ( verbosity > 0 )
			std::cout << "solved_tasks_count " << solved_tasks_count << std::endl;
		
		if ( process_sat_count ) {
			sat_count += process_sat_count;
			std::cout << "sat_count " << sat_count << std::endl;
			if ( finding_first_sat_time == 0 ) // first time only
				finding_first_sat_time = MPI_Wtime() - total_start_time;
		}
		
		WriteSolvingTimeInfo( solving_times, solved_tasks_count );
		
		if ( sat_count && !IsSolveAll )
			break; // exit if SAT set found
		
		if ( next_task_index < (int)all_tasks_count ) {
			// send new index of task
			MPI_Send( &next_task_index, 1, MPI_INT, current_status.MPI_SOURCE, 0, MPI_COMM_WORLD );
			// send to free core new task in format of minisat input masks
			copy( values_arr[next_task_index].begin(), values_arr[next_task_index].end(), mask_value );
			MPI_Send( mask_value, FULL_MASK_LEN, MPI_UNSIGNED, current_status.MPI_SOURCE, 0, MPI_COMM_WORLD );
			next_task_index++;
		}
	} // while ( solved_tasks_count < all_tasks_count )

	solving_iteration_count++;
	
	return true;
}

bool MPI_Solver :: ComputeProcessSolve()
{
	minisat22_wrapper m22_wrapper;
	Solver *S;
	MPI_Status status;
	std::stringstream sstream;
	int current_task_index;
	unsigned long long process_sat_count = 0;
	double cnf_time_from_node = 0.0;
	bool IsFirstTaskRecieved;
	std::string solving_info_file_name;
	std::string str;
	std::ifstream infile;
	int *var_choose_order_int;
	std::ifstream in;
	std::vector< std::vector<bool> > interrupted_problems_var_values_from_process, sat_assignment_from_process;
	char *char_send_array;
	int char_send_array_len = 0;
	
	for (;;) {
		if ( rank == 1 )
			std::cout << "new compute high level iteration" << std::endl;
		
		in.open( input_cnf_name ); // read every new batch because CNF can be changed
		Problem cnf;
		m22_wrapper.parse_DIMACS_to_problem(in, cnf);
		in.close();
		S = new Solver();
		S->addProblem(cnf);
		S->verbosity = 0;
		S->isPredict = false;
		
		IsFirstTaskRecieved = false;
		var_choose_order_int = new int[MAX_CORE_LEN];
		MPI_Recv( &core_len,                 1, MPI_INT,      0, 0, MPI_COMM_WORLD, &status );
		MPI_Recv( &all_tasks_count,          1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, &status );
		MPI_Recv( &solving_iteration_count,  1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status );
		MPI_Recv( &max_solving_time,         1, MPI_DOUBLE,   0, 0, MPI_COMM_WORLD, &status );
		MPI_Recv( &start_activity,           1, MPI_DOUBLE,   0, 0, MPI_COMM_WORLD, &status );
		MPI_Recv( var_choose_order_int,      MAX_CORE_LEN, MPI_INT, 0, 0, MPI_COMM_WORLD, &status );
		
		var_choose_order.resize(0);
		for( unsigned i=0; i < MAX_CORE_LEN; ++i ) {
			if ( var_choose_order_int[i] == -1 )
				break;
			if ( var_choose_order_int[i] <= 0 ) {
				std::cerr << "var_choose_order_int[i] <= 0" << std::endl;
				std::cerr << var_choose_order_int[i] << std::endl;
				return false;
			}
			var_choose_order.push_back( var_choose_order_int[i] );
		}
		delete[] var_choose_order_int;
		if ( rank == 1 ) {
			std::cout << "Received core_len "                << core_len                 << std::endl;
			std::cout << "Received all_tasks_count "         << all_tasks_count          << std::endl;
			std::cout << "Received solving_iteration_count " << solving_iteration_count  << std::endl;
			std::cout << "Received max_solving_time "        << max_solving_time         << std::endl;
			std::cout << "Received start_activity "          << start_activity           << std::endl;
			std::cout << "Received var_choose_order.size() " << var_choose_order.size()  << std::endl;
			for ( unsigned i=0; i < var_choose_order.size(); ++i )
				std::cout << var_choose_order[i] << " ";
			std::cout << std::endl;
		}
		
		S->core_len         = core_len;
		S->max_solving_time = max_solving_time;
		S->start_activity   = start_activity;
		S->resetVarActivity();
		
		if ( rank == 1 )
			std::cout << "before loop of recv tasks" << std::endl;
		for (;;) {
			// get index of current task
			MPI_Recv( &current_task_index, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status );
			
			if ( ( rank == 1 ) && ( verbosity > 2 ) )
				std::cout << "recv current_task_index " << current_task_index << std::endl;
			
			if ( current_task_index == -1 ) {
				if ( rank == 1 )
					std::cout << "breaking low level loop while computing" << std::endl;
				break; // stop and get new init values for solving iteration
			}
			else if ( current_task_index == -2 ) {
				MPI_Finalize( ); // finalize-message from control process
				break;
			}
			
			// with assumptions file we need only current_task_index for reading values from file 
			if ( !IsFirstTaskRecieved ) {
				MPI_Recv( full_mask, FULL_MASK_LEN, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, &status );
				MPI_Recv( part_mask, FULL_MASK_LEN, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, &status );
				IsFirstTaskRecieved = true;
			}
			MPI_Recv( mask_value, FULL_MASK_LEN, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, &status );
			
			if ( ( rank == 1 ) && ( verbosity > 2 ) ){
				std::cout << "full_mask" << std::endl;
				for( unsigned i=0; i<FULL_MASK_LEN; ++i )
					std::cout << full_mask[i] << " ";
				std::cout << std::endl;
				std::cout << "part_mask" << std::endl;
				for( unsigned i=0; i<FULL_MASK_LEN; ++i )
					std::cout << part_mask[i] << " ";
				std::cout << std::endl;
				std::cout << "mask_value" << std::endl;
				for( unsigned i=0; i<FULL_MASK_LEN; ++i )
					std::cout << mask_value[i] << " ";
				std::cout << std::endl;
			}
			
			if ( !SolverRun( S, process_sat_count, cnf_time_from_node, current_task_index, 
				             interrupted_problems_var_values_from_process, sat_assignment_from_process ) ) 
			{ 
				std::cout << std::endl << "Error in SolverRun"; 
				return false; 
			}
			
			if ( verbosity > 0 )
				std::cout << "process_sat_count is " << process_sat_count << std::endl;

			if ( ( verbosity > 0 ) && ( interrupted_problems_var_values_from_process.size() > 0 ) )
				std::cout << "interrupted_problems_var_values_from_process.size() " << interrupted_problems_var_values_from_process.size() << std::endl;
			
			MPI_Send( &process_sat_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD );
			MPI_Send( solving_times, SOLVING_TIME_LEN, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );

			if ( ( rank == 1 ) && ( verbosity > 0 ) ) {
				std::cout << "sended process_sat_count " << std::endl;
				std::cout << "sended solving_times " << std::endl;
			}
			
			// send info about interrupted problems
			if ( interrupted_problems_var_values_from_process.size() > 0 ) {
				char_send_array = new char[interrupted_problems_var_values_from_process.size() * interrupted_problems_var_values_from_process[0].size()];
				char_send_array_len = 0;
				for ( auto &x : interrupted_problems_var_values_from_process )
					for ( unsigned t = 0; t < x.size(); t++ )
						char_send_array[char_send_array_len++] = (x[t] == true ? '1' : '0');
			}
			else {
				char_send_array_len = 1;
				char_send_array = new char[char_send_array_len];
				char_send_array[0] = '2';
			}
			
			MPI_Send( char_send_array, char_send_array_len, MPI_CHAR, 0, 0, MPI_COMM_WORLD );
			delete[] char_send_array;
			
			// send satisfying assignments if such occur
			if ( sat_assignment_from_process.size() > 0 ) {
				char_send_array = new char[sat_assignment_from_process.size() * sat_assignment_from_process[0].size()];
				char_send_array_len = 0;
				for ( auto &x : sat_assignment_from_process )
					for ( unsigned t = 0; t < x.size(); t++ )
						char_send_array[char_send_array_len++] = (x[t] == true ? '1' : '0');
			}
			else {
				char_send_array_len = 1;
				char_send_array = new char[char_send_array_len];
				char_send_array[0] = '2';
			}

			MPI_Send( char_send_array, char_send_array_len, MPI_CHAR, 0, 0, MPI_COMM_WORLD );
			delete[] char_send_array;
			
			if ( ( rank == 1 ) && ( verbosity > 0 ) ) {
				std::cout << "char_send_arra " << std::endl;
				std::cout << "char_send_array_len " << char_send_array_len << std::endl;
			}
		}
		
		delete S;
	}
	
	return true;
}

//---------------------------------------------------------
bool MPI_Solver :: MPI_ConseqSolve( int argc, char **argv )
{
// read all cnf in current dir and start then in conseq mode
	std::string dir = std::string(".");
	int cnf_count = 0;
	std::vector<std::string> files;
	std::vector<std::string> cnf_files;
	std::fstream out_file;
	int current_obj_val = -1;
	unsigned long long process_sat_count= 0;
	double cnf_time_from_node;
	//PBSolver_cut pbs_cut;
	std::stringstream solve_sstream;
	double start_sec;
	double final_sec;
	Solver *S;
	std::vector<std::vector<bool>> interrupted_problems_var_values_from_process, sat_assignments_from_process;

	// MPI start
	//MPI_Request request;
	int corecount = 10, rank = 0;
	/*MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &corecount );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );*/

	files     = std::vector<std::string>( );
	cnf_files = std::vector<std::string>( );
	// get all files in current dir
	getdir( dir, files );

	for ( unsigned i = 0; i < files.size( ); i++ ) {
		if ( files[i].find( ".cnf" ) != std::string::npos ) {
			cnf_count++;
			cnf_files.push_back( files[i] );
			if ( rank == 0 ) 
				std::cout << std::endl << "founded cnf " << files[i].c_str( );
		}
	}

	if ( cnf_count > corecount ) {
		if ( rank == 0 ) {
			std::cout << std::endl << "Warning. Count of cnf-file > corecount";
			std::cout << std::endl << "Only first " << corecount << " cnf will be processed";
			std::cout << std::endl << "cnf_count changed to corecount";
		}
		cnf_count = corecount;
	}
	else if ( cnf_count == 0 ) {
		if ( rank == 0 ) std::cout << std::endl << "Error. No cnf-files in dir";
		return false;
	}

	if ( rank > cnf_count - 1 ) 
		std::cout << std::endl << "core # " << rank << " with no job";
	else {// do job
		std::stringstream sstream;
		sstream << "answer_" << rank + 1;
		std::string out_file_name = sstream.str( );
		sstream.str( "" );
		sstream.clear();
		
		start_sec = MPI_Wtime( ); // get init time
		input_cnf_name = &cnf_files[rank][0]; // set current name of file
		std::cout << std::endl << "input_cnf_name is " << input_cnf_name;
		unsigned int zero_mask[FULL_MASK_LEN];
		for ( int i = 0; i < FULL_MASK_LEN; i++ )
			zero_mask[i] = 0;
		if ( !ReadIntCNF( ) ) { // Read original CNF 
			std::cout << "\n Error in ReadIntCNF" << std::endl; return 1;
		}
		std::cout << std::endl << "end of ReadIntCNF";
		if ( rank == 0 ) 
			PrintParams( );
		if ( !IsPB ) {
			int current_task_index = 0;
			std::cout << std::endl << std::endl << "Standart mode of SAT solving";
			if ( !SolverRun( S, process_sat_count, cnf_time_from_node, current_task_index, 
				             interrupted_problems_var_values_from_process, sat_assignments_from_process ) ) 
			{
				std::cout << std::endl << "Error in SolverRun"; 
				return false;
			}
			if ( process_sat_count ) {
				if ( !AnalyzeSATset( cnf_time_from_node ) ) {
					// is't needed to deallocate memory - MPI_Abort will do it	
					std::cout << "\n Error in Analyzer" << std::endl;
					MPI_Abort( MPI_COMM_WORLD, 0 );
					return false;
				}
			}
		}
		
		final_sec = MPI_Wtime( ) - start_sec;
		
		sstream << input_cnf_name << " " << final_sec << " sec" << std::endl;
		out_file.open( out_file_name.c_str( ), std::ios_base :: out ); // open and clear out file
		out_file << sstream.rdbuf( );
		out_file << solve_sstream.rdbuf( );
		std::cout << std::endl << "*** sstream " << sstream.str( );
		std::cout << std::endl << "*** solve_sstream " << solve_sstream.str( );
		out_file.close( ); 		
	}

	MPI_Finalize( ); // MPI end

	std::cout << std::endl << "End of ControlConseqProcessSolve";
	return true;
}

void MPI_Solver :: PrintParams( )
{
	std::cout << std::endl << "solver_name is "           << solver_name;                         
	std::cout << std::endl << "koef_val is "              << koef_val;             		
	std::cout << std::endl << "schema_type is "		      << schema_type;           
	std::cout << std::endl << "full_mask_var_count is "   << full_mask_var_count;  
	std::cout << std::endl << "proc_count is "            << corecount;            
	std::cout << std::endl << "core_len is "              << core_len;                 
	std::cout << std::endl << "start_activity is "        << start_activity;
	std::cout << std::endl << "IsConseq is "              << IsConseq;         
	std::cout << std::endl << "IsPB is "                  << IsPB;                     
	std::cout << std::endl << "best_lower_bound is "      << best_lower_bound;        
	std::cout << std::endl << "upper_bound is "	          << upper_bound;			  
	std::cout << std::endl << "PB_mode is "               << PB_mode;
	std::cout << std::endl << "exch_activ is "            << exch_activ;
	std::cout << std::endl << "IsSolveAll exch_activ is " << IsSolveAll;
	std::cout << std::endl << "max_solving_time "         << max_solving_time;
	std::cout << std::endl << "max_nof_restarts "         << max_nof_restarts;
	std::cout << std::endl << "max_solving_time_koef "    << max_solving_time_koef;
	std::cout << std::endl << "var_count "                << var_count;
	std::cout << std::endl;
}

//---------------------------------------------------------
bool MPI_Solver :: WriteTimeToFile( double time_sec )
{
// Write time of solving to output file
	std::string line_buffer,
		   TimeIs_char = "\n Time: ",
		   new_line_char = "\n",
		   hour_char = " h ",
		   minute_char = " m ",
		   second_char = " s ";
	int real_hours = -1,
		real_minutes = -1,
		real_seconds = -1;

	std::cout << std::endl << "That took %f seconds " << time_sec << std::endl;
	
	cpuTimeInHours( time_sec, real_hours, real_minutes, real_seconds );
	std::stringstream sstream;
	sstream << TimeIs_char << real_hours << hour_char 
		    << real_minutes << minute_char << real_seconds << second_char;
	
	std::ofstream solving_info_file;
	solving_info_file.open( solving_info_file_name.c_str( ), std::ios::app );
	if ( !solving_info_file.is_open( ) ) {
		std :: cout << "Error in opening of output file " << solving_info_file_name << std::endl;
		return false;
	}
	solving_info_file << sstream.rdbuf( );
	solving_info_file.close( );

	return true;
}