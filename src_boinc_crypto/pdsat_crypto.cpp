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
string previous_results_str = "";

int bivium_template_array[] = {
#include "bivium_template.inc"
};

bool do_work( string &input_path, string &final_result_str );
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
    /*ifstream infile( input_path.c_str() );
    if ( !infile.is_open() ) {
		fprintf(stderr, "%s APP: app infile open failed:\n",
            boinc_msg_prefix(buf, sizeof(buf))
        );
        exit(-1);
    }*/
	
	// See if there's a valid checkpoint file.
    boinc_resolve_filename_s( CHECKPOINT_FILE, chpt_path );
	chpt_file.open( chpt_path );
	if ( !chpt_file.is_open() ) {
		fprintf(stderr, "%s APP: app chpt open failed:\n",
            boinc_msg_prefix(buf, sizeof(buf))
        );
        exit(-1);
    }
	
	string str;
	chpt_file >> last_iteration_done >> total_problems_count;
	while ( getline( chpt_file, str ) )
		previous_results_str += str;
    chpt_file.close();

	string final_result_str;
	if ( !do_work( input_path, final_result_str ) ) {
		fprintf( stderr, "%s APP: do_work() failed:\n" );
		perror("do_work");
        exit(1);
	}

	// resolve and open output file
    boinc_resolve_filename_s( OUTPUT_FILENAME, output_path );
	outfile.open( output_path.c_str() );
    if ( !outfile.is_open() ) {
        fprintf(stderr, "%s APP: app output open failed:\n",
            boinc_msg_prefix(buf, sizeof(buf))
        );
		exit(-1);
    }
	
	outfile << final_result_str;
	outfile.close();

    boinc_finish(0);
}

bool do_work( string &input_path, string &final_result_str )
{
	int retval;
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
	fprintf( stderr, " %d ", S->max_nof_restarts );
	S->addProblem( cnf ); 

	// read assignments from input file
	MPI_Base mpi_b;
	vec< vec<Lit> > dummy_vec;
	int current_task_index = 0;
	mpi_b.known_assumptions_file_name = input_path;
	mpi_b.MakeAssignsFromFile( current_task_index, dummy_vec );

	double current_time = 0;
	double time_last_checkpoint = Minisat :: cpuTime();

	string current_result_str = "";
	if ( previous_results_str.find(" SAT") != string :: npos ) {
		current_result_str = previous_results_str;
	}

		// solve

		current_time = Minisat :: cpuTime() - time_last_checkpoint;
		// skip some time in case of fast checkpoint to make it correct
		if ( current_time >= time_last_checkpoint + MIN_CHECKPOINT_INTERVAL_SEC ) {
			//current_result_str = previous_results_str + ;
			// checkpoint current position and results
			//if ( ( boinc_is_standalone() ) || ( boinc_time_to_checkpoint() ) ) {
				//retval = do_checkpoint( ls.problems_solved + last_iteration_done, total_problems_count, current_result_str );
				if (retval) {
					fprintf(stderr, "APP: checkpoint failed %d\n", retval );
					exit(retval);
				}
				boinc_checkpoint_completed();
			//}
			time_last_checkpoint = current_time;
		}

	delete S;

	final_result_str = current_result_str;
	
	return true;
}

int do_checkpoint( unsigned current_solved, unsigned total_tasks, string &current_result_str ) {
    int retval;
    string resolved_name;

    FILE* f = fopen("temp", "w");
    if (!f) return 1;
    fprintf( f, "%d %d \n%s", current_solved, total_tasks, current_result_str.c_str() );
    fclose( f );

    boinc_resolve_filename_s(CHECKPOINT_FILE, resolved_name);
    retval = boinc_rename( "temp", resolved_name.c_str() );
    if ( retval ) return retval;

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