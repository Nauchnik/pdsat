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

unsigned long long last_iteration_done = 0;
unsigned long long total_problems_count = 0;
string previous_results_str = "";

int diag10_2_cnf_array[] = {
#include "diag10_2.inc"
};

bool do_work( string &input_path, string &current_result_str );
int do_checkpoint( unsigned long long current_solved, unsigned long long total_tasks, string &current_result_str );

int main(int argc, char **argv) {
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
		chpt_file >> last_iteration_done >> total_problems_count;
		while ( getline( chpt_file, str ) )
			previous_results_str += str;
		chpt_file.close();
	}
	
	string current_result_str = ""; // if SAT then it will be changed
	if ( !do_work( input_path, current_result_str ) ) {
		fprintf( stderr, "APP: do_work() failed:\n" );
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
	latin_square ls;
	int retval;
	string error_msg;
	
	fprintf(stderr, "input_path: %s\n", input_path.c_str());

	if ( !ls.ReadLiteralsFromFile( input_path, error_msg ) ) {
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

	vector< vector<int> > :: iterator positive_literals_it;
	double current_time = 0, time_last_checkpoint = 0;
	clock_t clk_start = clock( );
	time_last_checkpoint = (double)(clock( ) - clk_start)/(double)(CLOCKS_PER_SEC);
	
	fprintf(stderr, " restarts %d \n", MAX_NOF_RESTARTS );
	fprintf(stderr, " ls.positive_literals size %d \n", ls.positive_literals.size() );
	
	for ( positive_literals_it = ls.positive_literals.begin() + last_iteration_done; 
		  positive_literals_it != ls.positive_literals.end(); positive_literals_it++ ) 
	{
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
		/*while ( current_time < time_last_checkpoint + MIN_CHECKPOINT_INTERVAL_SEC ) 
			current_time = (double)(clock( ) - clk_start)/(double)(CLOCKS_PER_SEC); 
		time_last_checkpoint = current_time;*/
		current_result_str = previous_results_str + ls.all_answers.str();
		
		// checkpoint current position and results
		//if ( ( boinc_is_standalone() ) || ( boinc_time_to_checkpoint() ) ) {
			retval = do_checkpoint( (unsigned long long)(ls.problems_solved) + last_iteration_done, total_problems_count, current_result_str );
            if (retval) {
                fprintf(stderr, "APP: checkpoint failed %d\n", retval );
                exit(retval);
            }
			boinc_checkpoint_completed();
        //}
		
		delete S;
	}
	
	return true;
}

int do_checkpoint( unsigned long long current_solved, unsigned long long total_tasks, string &current_result_str ) 
{
	int retval;
    string resolved_name;
	
	ofstream temp_ofile( "temp" );
	if ( !temp_ofile.is_open() ) 
		return 1;
    temp_ofile << current_solved << " " << total_tasks << " " << current_result_str;
    temp_ofile.close();
	
    boinc_resolve_filename_s(CHECKPOINT_FILE, resolved_name);
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