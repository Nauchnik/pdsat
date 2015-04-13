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

#include "latin_squares.h"

using namespace std;

#define CHECKPOINT_FILE "chpt"
#define INPUT_FILENAME "in"
#define OUTPUT_FILENAME "out"

const int MAX_NOF_RESTARTS            = 5000;
const int MIN_CHECKPOINT_INTERVAL_SEC = 2;
const int POSITIVE_LITERALS_MIN_SIZE  = 8;

int last_iteration_done = 0;
int total_problems_count = 0;
string previous_results_str = "";

int diag10_2_cnf_array[] = {
#include "diag10_2.inc"
};

bool do_work( FILE *infile, string &final_result_str );
int do_checkpoint( unsigned current_solved, unsigned total_tasks, string &final_result_str );

int main(int argc, char **argv) {
    int nchars = 0, retval, n;
    double fsize;
    char input_path[512], output_path[512], chkpt_path[512], buf[256];
    MFILE outfile;
    FILE *chpt_file, *infile;
	
    retval = boinc_init();
    if (retval) {
        fprintf(stderr, "%s boinc_init returned %d\n",
            boinc_msg_prefix(buf, sizeof(buf)), retval
        );
        exit(retval);
    }

    // open the input file (resolve logical name first)
    boinc_resolve_filename(INPUT_FILENAME, input_path, sizeof(input_path));
    infile = boinc_fopen(input_path, "r");
    if (!infile) {
        fprintf(stderr,
            "%s Couldn't find input file, resolved name %s.\n",
            boinc_msg_prefix(buf, sizeof(buf)), input_path
        );
        exit(-1);
    }

    // get size of input file (used to compute fraction done)
    file_size( input_path, fsize );

	// resolve output file
    boinc_resolve_filename(OUTPUT_FILENAME, output_path, sizeof(output_path));
	
	// See if there's a valid checkpoint file.
	char string_input[1024];
    boinc_resolve_filename(CHECKPOINT_FILE, chkpt_path, sizeof(chkpt_path));
    chpt_file = boinc_fopen(chkpt_path, "r");
	
    if ( chpt_file ) {
        n = fscanf( chpt_file, "%d %d", &last_iteration_done, &total_problems_count );
		while ( fgets(string_input, 1024, chpt_file ) )
			previous_results_str += string_input;
        fclose( chpt_file );
    }
	
	retval = outfile.open(output_path, "w");
    if (retval) {
        fprintf(stderr, "%s APP: app output open failed:\n",
            boinc_msg_prefix(buf, sizeof(buf))
        );
        fprintf(stderr, "%s resolved name %s, retval %d\n",
            boinc_msg_prefix(buf, sizeof(buf)), output_path, retval
        );
        perror("open");
        exit(1);
    }

	string final_result_str;
	if ( !do_work( infile, final_result_str ) ) {
		fprintf(stderr, "%s APP: do_work() failed:\n" );
		perror("do_work");
        exit(1);
	}

	outfile.puts( final_result_str.c_str() );
    retval = outfile.flush();
    if (retval) {
        fprintf( stderr, "%s APP: failed %d\n",
                 boinc_msg_prefix(buf, sizeof(buf)), retval
        );
        exit(1);
    }

    boinc_finish(0);
}

bool do_work( FILE *infile, string &final_result_str )
{
	latin_square ls;
	int retval;
	string error_msg;
	
	if ( !ls.ReadLiteralsFromFile( infile, error_msg ) ) {
		fprintf(stderr, "APP: error in ls.ReadLiteralsFromFile\n");
		fprintf(stderr, "APP: %s\n", error_msg.c_str());
		return false;
	}

	if ( ls.positive_literals.size() == 0 ) {
		fprintf(stderr, "APP: ls.positive_literals.size() == 0\n");
		return false;
	}

	if ( ls.positive_literals[0].size() < POSITIVE_LITERALS_MIN_SIZE ) {
		fprintf(stderr, "APP: ls.positive_literals[0].size() < POSITIVE_LITERALS_MIN_SIZE\n");
		fprintf(stderr, "%d < %d\n", ls.positive_literals[0].size(), POSITIVE_LITERALS_MIN_SIZE);
		return false;
	}

	if ( total_problems_count == 0 ) // if first launch and no chkpt file
		total_problems_count = ls.positive_literals.size();
	
	if ( ls.problem_type == "diag" ) {
		ls.cnf_array.resize( sizeof(diag10_2_cnf_array)  / sizeof(diag10_2_cnf_array[0]) );
		for ( unsigned i = 0; i < ls.cnf_array.size(); i++ ) 
			ls.cnf_array[i] = diag10_2_cnf_array[i];
		fprintf(stderr, " ls.problem_type diag\n ");
	}
	else {
		fprintf(stderr, " not diag ");
		return false;
	}
	
	minisat22_wrapper m22_wrapper;
	Problem cnf;
	m22_wrapper.parse_DIMACS_from_inc( ls.cnf_array, cnf );
	Solver *S;

	string current_result_str = "";
	vector< vector<int> > :: iterator positive_literals_it;
	double current_time = 0, time_last_checkpoint = 0;
	clock_t clk_start = clock( );
	time_last_checkpoint = (double)(clock( ) - clk_start)/(double)(CLOCKS_PER_SEC);
	
	fprintf(stderr, " restarts %d \n", MAX_NOF_RESTARTS );
	fprintf(stderr, " ls.positive_literals size %d \n", ls.positive_literals.size() );
	
	for ( positive_literals_it = ls.positive_literals.begin() + last_iteration_done; 
		  positive_literals_it != ls.positive_literals.end(); positive_literals_it++ ) 
	{
		//for ( unsigned i = 0; i < positive_literals_it->size(); i++ )
		//	fprintf(stderr, " %d ", positive_literals_it->at(i) );
		//fprintf(stderr, " \n " );

		S = new Solver();
		S->max_nof_restarts = MAX_NOF_RESTARTS;
		S->problem_type = ls.problem_type;

		S->addProblem( cnf ); // add initial CNF every time
		if ( !ls.SolveOneProblem( S, positive_literals_it, clk_start ) ) {
			fprintf(stderr, "Interrupted in client.cpp after SolveOneProblem() returned false");
			break;
		}
		current_time = (double)(clock( ) - clk_start)/(double)(CLOCKS_PER_SEC);
		// skip some time in case of fast checkpoint to make it correct
		while ( current_time < time_last_checkpoint + MIN_CHECKPOINT_INTERVAL_SEC ) 
			current_time = (double)(clock( ) - clk_start)/(double)(CLOCKS_PER_SEC); 
		time_last_checkpoint = current_time;
		current_result_str = previous_results_str + ls.all_answers.str();

		// checkpoint current position and results
		//if ( ( boinc_is_standalone() ) || ( boinc_time_to_checkpoint() ) ) {
			retval = do_checkpoint( ls.problems_solved + last_iteration_done, total_problems_count, current_result_str );
            if (retval) {
                fprintf(stderr, "APP: checkpoint failed %d\n", retval );
                exit(retval);
            }
			boinc_checkpoint_completed();
        //}

		delete S;
	}
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