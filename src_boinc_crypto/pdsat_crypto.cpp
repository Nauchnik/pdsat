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

using namespace std;

#define CHECKPOINT_FILE "chpt"
#define INPUT_FILENAME  "in"
#define OUTPUT_FILENAME "out"

const int MAX_NOF_RESTARTS            = 5000;
const int MIN_CHECKPOINT_INTERVAL_SEC = 10;

int last_iteration_done = 0;
int total_problems_count = 0;
string previous_results_str;

int bivium_template_array[] = {
#include "bivium_template.inc"
};

bool do_work( string &input_path, string &current_result_str );
int do_checkpoint( unsigned current_solved, unsigned total_tasks, string &final_result_str );

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
    ofstream outfile;
    fstream chpt_file;
	
    // open the input file (resolve logical name first)
	boinc_resolve_filename_s( INPUT_FILENAME, input_path );
	
	// See if there's a valid checkpoint file.
    boinc_resolve_filename_s( CHECKPOINT_FILE, chpt_path );
	chpt_file.open( chpt_path, ios_base :: out );
	if ( !chpt_file.is_open() ) {
		fprintf(stderr, "%s APP: app chpt open failed:\n",
            boinc_msg_prefix(buf, sizeof(buf))
        );
        exit(-1);
    }
	string str;
	chpt_file >> last_iteration_done >> total_problems_count;
	previous_results_str = "";
	while ( getline( chpt_file, str ) )
		previous_results_str += str;
    chpt_file.close();

	string current_result_str; // if SAT then it will be changed
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
	string error_msg;
	string problem_type;
	// before assignments there are option string strating with "h"
	ifstream infile( input_path.c_str() );
	infile >> problem_type;
	vector<int> cnf_array;
	infile.close();
	
	fprintf( stderr, problem_type.c_str() );
	if ( problem_type.find( "bivium" ) != string::npos ) {
		cnf_array.resize( sizeof(bivium_template_array)  / sizeof(bivium_template_array[0]) );
		for ( unsigned i = 0; i < cnf_array.size(); ++i ) 
			cnf_array[i] = bivium_template_array[i];
	}

	// read initial CNF from structure and add it to Solver
	minisat22_wrapper m22_wrapper;
	Problem cnf;
	m22_wrapper.parse_DIMACS_from_inc( cnf_array, cnf );
	Solver *S = new Solver();
	S->max_nof_restarts = MAX_NOF_RESTARTS;
	S->verbosity = 0;
	fprintf( stderr, " %d ", S->max_nof_restarts );
	S->addProblem( cnf ); 

	// read assignments from input file
	MPI_Base mpi_b;
	ifstream ifile( input_path );
	string str, str1;
	stringstream sstream;
	int val;
	while ( getline( ifile, str ) ) {
		sstream << str;
		sstream >> str1;
		if ( str1 == "var_set" ) {
			while ( sstream >> val )
				mpi_b.var_choose_order.push_back( val );
			sort( mpi_b.var_choose_order.begin(), mpi_b.var_choose_order.end() );
			break;
		}
		sstream.clear(); sstream.str();
	}
	ifile.close();
	vec< vec<Lit> > dummy_vec;
	int current_task_index = 0;
	mpi_b.known_assumptions_file_name = input_path;
	mpi_b.MakeAssignsFromFile( current_task_index, dummy_vec );
	
	double time_last_checkpoint = Minisat :: cpuTime();
	double current_time = 0;

	if ( previous_results_str.find(" SAT") != string :: npos )
		current_result_str = previous_results_str;
	else {
		current_result_str = problem_type + " var_set size " + inttostr(mpi_b.var_choose_order.size()) + ": ";
		for ( unsigned i=0; i < mpi_b.var_choose_order.size(); ++i )
			current_result_str += (inttostr( mpi_b.var_choose_order[i] ) + " ");
	}
	
	lbool ret;
	int retval;
	int current_launch_problems_solved = 0;
	total_problems_count = dummy_vec.size();
	double one_solving_time;
	bool isSAT = false;
	
	for ( int i = last_iteration_done; i < dummy_vec.size(); ++i ) {
		S->last_time = Minisat :: cpuTime();
		ret = S->solveLimited( dummy_vec[i] );
		one_solving_time = Minisat :: cpuTime() - S->last_time;
		current_launch_problems_solved++;
		
		if ( ret == l_True ) {
			current_result_str += " SAT ";
			for ( int i=0; i < S->model.size(); i++ )
				current_result_str += ( S->model[i] == l_True) ? '1' : '0';
			isSAT = true;
		}
		
		current_time = Minisat :: cpuTime() - time_last_checkpoint;
		// skip some time in case of fast checkpoint to make it correct
		if ( current_time >= time_last_checkpoint + MIN_CHECKPOINT_INTERVAL_SEC ) {
			// checkpoint current position and results
			//if ( ( boinc_is_standalone() ) || ( boinc_time_to_checkpoint() ) ) {
				retval = do_checkpoint( current_launch_problems_solved + last_iteration_done, total_problems_count, current_result_str );
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
	
	if ( !isSAT )
		current_result_str += "UNSAT";
	
	return true;
}

int do_checkpoint( unsigned current_solved, unsigned total_tasks, string &current_result_str ) {
    int retval;
    string resolved_name;

    ofstream temp_ofile( "temp" );
	if ( !temp_ofile.is_open() ) 
		return 1;
    temp_ofile << current_solved << " " << total_tasks << " " << current_result_str;
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