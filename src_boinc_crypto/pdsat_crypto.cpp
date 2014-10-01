// BOINC client application PDSAT

#ifdef _WIN32
#include "boinc_win.h"
#else
#include "config.h"
#include <cstdio>
#include <cctype>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <csignal>
#include <unistd.h>
#endif

#include "str_util.h"
#include "util.h"
#include "filesys.h"
#include "boinc_api.h"
#include "mfile.h"

#include <fstream>
#include "minisat22_wrapper.h"
#include "mpi_base.h"
#include "dminisat/dminisat.h"

using namespace std;

#define CHECKPOINT_FILE "chpt"
#define INPUT_FILENAME  "in"
#define OUTPUT_FILENAME "out"

const int MAX_NOF_RESTARTS            = 5000;
const double MAX_SOLVING_TIME         = 360.0;
const int MIN_CHECKPOINT_INTERVAL_SEC = 10;

unsigned long long last_iteration_done = 0;
unsigned long long total_problems_count = 0;
string previous_results_str;
double max_solved_time = 0.0;

int bivium_template_array[] = {
#include "bivium_template.inc"
};

int a5_1_114_template_array[] = {
#include "a5_1_114_template.inc"
};

bool do_work( string &input_path, string &current_result_str );
int do_checkpoint( unsigned long long current_solved, unsigned long long total_tasks, double max_solved, string &final_result_str );

int main( int argc, char **argv ) {
    char buf[256];
	int retval = boinc_init();
    if ( retval ) {
        fprintf(stderr, "%s boinc_init returned %d\n",
            boinc_msg_prefix(buf, sizeof(buf)), retval
        );
        exit( retval );
    }

	string input_path, output_path, chpt_path;
	string str;
    ofstream outfile;
    fstream chpt_file;
	
    // open the input file (resolve logical name first)
	boinc_resolve_filename_s( INPUT_FILENAME, input_path );
	
	// See if there's a valid checkpoint file.
	previous_results_str = "";
    boinc_resolve_filename_s( CHECKPOINT_FILE, chpt_path );
	chpt_file.open( chpt_path.c_str(), ios_base :: in );
	if ( chpt_file.is_open() ) {
		chpt_file >> last_iteration_done >> total_problems_count >> max_solved_time;
		while ( getline( chpt_file, str ) )
			previous_results_str += str;
		chpt_file.close();
	}
	
	string current_result_str = ""; // if SAT then it will be changed
	if ( !do_work( input_path, current_result_str ) ) {
		fprintf( stderr, "%s APP: do_work() failed:\n" );
		perror("do_work");
        exit(1);
	}

	// resolve, open and write answer to output file
    boinc_resolve_filename_s( OUTPUT_FILENAME, output_path );
	outfile.open( output_path.c_str() );
    if ( !outfile.is_open() ) {
        fprintf(stderr, "%s APP: app output open failed:\n",
            boinc_msg_prefix(buf, sizeof(buf))
        );
		exit(-1);
    }
	outfile << current_result_str;
	outfile.close();
	
    boinc_finish(0);
}

bool do_work( string &input_path, string &current_result_str )
{
	MPI_Base mpi_b;
	stringstream sstream;
	string str_tmp;
	bool isRangeMode = false;
	unsigned long long range1 = 0, range2 = 0;
	
	// read var_choose_order and assignments from file in text mode
	string problem_type, str, word1;
	int val;
	vector<int> known_var_values_vec;
	ifstream ifile( input_path.c_str() );
	getline( ifile, problem_type );
	fprintf( stderr, problem_type.c_str() );

	stringstream range_sstream;

	while ( getline( ifile, str ) ) {
		if ( str.find( "before_assignments") != string::npos )
			break;
		if ( str.find( "before_range") != string::npos ) {
			isRangeMode = true;
			getline( ifile, str );
			sstream << str;
			sstream >> range1 >> range2;
			sstream.str(""); sstream.clear();
			if ( range2 < range1 ) {
				cerr << "range2 < range1" << endl;
				return false;
			}
			range_sstream << "range1 " << range1 << " range2 " << range2;
			fprintf( stderr, " range1 %llu range2 %llu ", range1, range2 );
			break;
		}
		sstream << str;
		sstream >> word1;
		if ( word1 == "var_set" ) {
			while ( sstream >> val )
				mpi_b.var_choose_order.push_back( val );
			sort( mpi_b.var_choose_order.begin(), mpi_b.var_choose_order.end() );
		}
		else if ( isNumberOrMinus( word1[0] ) )
			known_var_values_vec.push_back( strtoint( word1 ) );
		sstream.clear(); sstream.str("");
	}
	ifile.close();
	
	fprintf( stderr, " mpi_b.var_choose_order.size() %d", mpi_b.var_choose_order.size() );
	fprintf( stderr, " var_values_vec.size() %d", known_var_values_vec.size() );
	
	if ( known_var_values_vec.size() == 0 ) {
		cerr << "var_values_vec.size == 0" << endl;
		return false;
	}
	
	// read initial CNF from structure and add it to Solver
	vector<int> cnf_array;
	if ( problem_type.find( "bivium" ) != string::npos ) {
		mpi_b.input_var_num = 177;
		cnf_array.resize( sizeof(bivium_template_array) / sizeof(bivium_template_array[0]) );
		for ( unsigned i = 0; i < cnf_array.size(); ++i ) 
			cnf_array[i] = bivium_template_array[i];
		fprintf( stderr, " bivium_template " );
	}
	else if ( problem_type.find( "a5_1" ) != string::npos ) {
		mpi_b.input_var_num = 64;
		cnf_array.resize( sizeof(a5_1_114_template_array) / sizeof(a5_1_114_template_array[0]) );
		for ( unsigned i = 0; i < cnf_array.size(); ++i ) 
			cnf_array[i] = a5_1_114_template_array[i];
		fprintf( stderr, " a5_1_114_template " );
	}
	else {
		fprintf( stderr, " problem_type != bivium, problem_type != a5_1" );
		exit(1);
	}
	minisat22_wrapper m22_wrapper;
	Problem cnf;
	m22_wrapper.parse_DIMACS_from_inc( cnf_array, cnf );
	Solver *S = new Solver();
	S->addProblem( cnf );
	S->max_nof_restarts = MAX_NOF_RESTARTS;
	S->max_solving_time = MAX_SOLVING_TIME;
	fprintf( stderr, " S->max_nof_restarts %d ", S->max_nof_restarts );
	fprintf( stderr, " S->max_solving_time %f ", S->max_solving_time );
	S->verbosity = 0;
	//S->cur_hack_type = hack_minigolf;
	// minisat core mod
	S->core_len = mpi_b.input_var_num;
	S->start_activity = 1;
	S->resetVarActivity();
	fprintf( stderr, " S->core_len %d ", S->core_len );
	fprintf( stderr, " S->start_activity %f ", S->start_activity );
	for( unsigned i = 0; i < cnf.size(); ++i )
		delete cnf[i];
	cnf.clear();
	//vec< vec<Lit> > dummy_vec;

	/*if ( !isRangeMode  ) {
		// find size of text block before bynary block
		ifile.open( input_path.c_str(), ios_base :: in | ios_base :: binary );
		ifile.seekg (0, ifile.end);
		int length = ifile.tellg();
		if ( length <= 0 ) {
			fprintf( stderr, " length of buffer is %d ", length );
			return false;
		}
		ifile.seekg (0, ifile.beg);
		char *buffer = new char[length + 1];
		buffer[length] = '\0';
		ifile.read( buffer, length ); // read text data as a block:
		char *result = strstr( buffer, "before_assignments" );
		int before_binary_length = -1;
		if ( result ) {
			while ( !isNumber( result[0] ) )
				result++;
			before_binary_length = result - buffer;
		}
		delete[] buffer;
		if ( before_binary_length <= 0 ) {
			fprintf( stderr, "before_binary_length <= 0");
			return false;
		}
		ifile.close();
		fprintf( stderr, " size of before_binary_length %d", before_binary_length );
		// make vector of assunptions basing on bunary data
		ifile.open( input_path.c_str(), ios_base :: in | ios_base :: binary );
		buffer = new char[before_binary_length + 1];
		buffer[before_binary_length] = '\0';
		ifile.read( buffer, before_binary_length );
		delete[] buffer;
		short int si;
		unsigned long long ul;
		ifile.read( (char*)&si, sizeof(si) ); // read header
		fprintf( stderr, " binary prefix %d", si );
		mpi_b.assumptions_count = 0;
		//while ( ifile.read( (char*)&ul, sizeof(ul) ) )
		//	mpi_b.assumptions_count++;
		ifile.clear();
		ifile.seekg( 0, ifile.end );
		long long int total_byte_length = ifile.tellg();
		ifile.close();
		mpi_b.assumptions_count = (total_byte_length - before_binary_length - 2) / sizeof(ul);
		fprintf( stderr, " mpi_b.assumptions_count %d", mpi_b.assumptions_count );
		mpi_b.known_assumptions_file_name = input_path;
		mpi_b.all_tasks_count = 1;
		int current_task_index = 0;
		if ( !mpi_b.MakeAssignsFromFile( current_task_index, before_binary_length, dummy_vec ) ) {
			cerr << "MakeAssignsFromFile()";
			return false;
		}
		fprintf( stderr, " vector of assumptions was made " );
	}*/
	
	double time_last_checkpoint = Minisat :: cpuTime();
	double current_time = 0;

	if ( previous_results_str.find(" SAT") != string :: npos )
		current_result_str = previous_results_str;
	else {
		sstream << mpi_b.var_choose_order.size();
		current_result_str = problem_type + " " + range_sstream.str() + " var_set size " + sstream.str() + ": ";
		sstream.clear(); sstream.str("");
		for ( unsigned i=0; i < mpi_b.var_choose_order.size(); ++i ) {
			sstream << mpi_b.var_choose_order[i];
			current_result_str += (sstream.str() + " ");
			sstream.clear(); sstream.str("");
		}
	}
	
	lbool ret;
	int retval;
	unsigned long long current_launch_problems_solved = 0;
	total_problems_count = range2-range1 + 1;
	double one_solving_time;
	bool isSAT = false;
	
	fprintf( stderr, "start of range mode " );
	fprintf( stderr, " before loop of solving " );
	double total_solving_time = Minisat :: cpuTime();
	
	//string core_activity_str;
	//sstream << "core_activity " << endl;

	int cur_var_ind;
	// read range of values and made assumptions

	boost::dynamic_bitset<> d_b;
	d_b.resize( mpi_b.var_choose_order.size() );
	vec<Lit> dummy;
	unsigned dummy_size = mpi_b.var_choose_order.size() + known_var_values_vec.size();
	dummy.resize( dummy_size );
	unsigned dummy_index;

	for ( unsigned long long i = last_iteration_done; i < total_problems_count; ++i ) {
		// make vector of assignments
		UllongToBitset( i, d_b );
		if ( d_b.size() > mpi_b.var_choose_order.size() ) {
			fprintf( stderr, "d_b.size() > mpi_b.var_choose_order.size()" );
			return false;
		}
		dummy_index = 0;
		for ( unsigned j=0; j < d_b.size(); ++j ) {
			cur_var_ind = mpi_b.var_choose_order[j] - 1;
			dummy[dummy_index++] = ( d_b[j] == 1 ) ? mkLit( cur_var_ind ) : ~mkLit( cur_var_ind );
		}
		// add known data to assumptions (initially in oneliteral clauses)
		for ( unsigned j=0; j < known_var_values_vec.size(); ++j ) {
			cur_var_ind = abs( known_var_values_vec[j] ) - 1;
			dummy[dummy_index++] = known_var_values_vec[j] > 0 ? mkLit( cur_var_ind ) : ~mkLit( cur_var_ind );
		}
		
		S->last_time = Minisat :: cpuTime();
		ret = S->solveLimited( dummy );
		//fprintf( stderr, "\n S->solveLimited() done " );
		one_solving_time = Minisat :: cpuTime() - S->last_time;
		if ( max_solved_time < one_solving_time )
			max_solved_time = one_solving_time;
		current_launch_problems_solved++;
		
		if ( ret == l_True ) {
			current_result_str += " SAT ";
			for ( int i=0; i < S->model.size(); i++ )
				current_result_str += ( S->model[i] == l_True) ? '1' : '0';
			current_result_str += " ";
			isSAT = true;
		}
		
		current_time = Minisat :: cpuTime();
		// skip some time in case of fast checkpoint to make it correct
		if ( current_time >= time_last_checkpoint + MIN_CHECKPOINT_INTERVAL_SEC ) {
			// checkpoint current position and results
			//if ( ( boinc_is_standalone() ) || ( boinc_time_to_checkpoint() ) ) {
				retval = do_checkpoint( current_launch_problems_solved + last_iteration_done, total_problems_count, max_solved_time, current_result_str );
				if ( retval ) {
					fprintf(stderr, "APP: checkpoint failed %d\n", retval );
					exit( retval );
				}
				boinc_checkpoint_completed();
			//}
			time_last_checkpoint = Minisat :: cpuTime();
		}
	}
	delete S;

	total_solving_time = Minisat :: cpuTime() - total_solving_time;
	
	unsigned long long total_solved = current_launch_problems_solved + last_iteration_done;
	sstream << total_solved;
	current_result_str += "solved " + sstream.str() + " ";
	sstream.clear(); sstream.str("");
	sstream << total_solving_time;
	current_result_str += "total_time " + sstream.str() + " ";
	sstream.clear(); sstream.str("");
	sstream << max_solved_time;
	current_result_str += "max_time " + sstream.str() + " ";
	sstream.clear(); sstream.str("");
	if ( !isSAT )
		current_result_str += "UNSAT";
	
	return true;
}

int do_checkpoint( unsigned long long current_solved, unsigned long long total_tasks, 
				   double max_solving_time, string &current_result_str ) {
    int retval;
    string resolved_name;
	
    ofstream temp_ofile( "temp" );
	if ( !temp_ofile.is_open() ) 
		return 1;
    temp_ofile << current_solved << " " << total_tasks << " " << max_solved_time << " " << current_result_str;
    temp_ofile.close();

    boinc_resolve_filename_s( CHECKPOINT_FILE, resolved_name );
    retval = boinc_rename( "temp", resolved_name.c_str() );
    if ( retval ) 
		return retval;

	boinc_fraction_done( ( double )current_solved / ( double )total_tasks );

    return 0;
}

#ifdef _WIN32
int WINAPI WinMain(HINSTANCE hInst, HINSTANCE hPrevInst, LPSTR Args, int WinMode) {
    LPSTR command_line;
    char* argv[100];
    int argc;

    command_line = GetCommandLine();
    argc = parse_command_line( command_line, argv );
    return main(argc, argv);
}
#endif