#ifndef _WIN32
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#endif

#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <ctype.h>
#include <fstream>
#include <string>
#include <sstream>

#include "latin_squares.h"

const int MAX_NOF_RESTARTS = 1000;
const int PROBLEMS_IN_WU = 20;

#pragma warning( disable : 4996 )

#define INPUT_FILENAME	"in"
#define OUTPUT_FILENAME	"out"

int diag10_2_cnf_array[] = {
#include "../src_common/ls10_2_diag.inc"
};

fstream infile;
fstream outfile;
ifstream chptfile;
string infile_name;
string outfile_name;
string ckptfile_name;

int main( int argc, char* argv[] )
{
	string str;
	stringstream sstream;
	latin_square ls;
	ifstream tasks_file, values_file;
	vector<int> tasks_vec;
	vector< vector<int> > known_values_vec;
	bool IsTasksFile = false;
	bool IsValuesFile = false;
	int int_val;
	
#ifdef _DEBUG
	argc = 5;
	argv[1] = "10";
	argv[2] = "2";
	argv[3] = "5";
	argv[4] = "diag";
	//argv[5] = "values.txt";
	//argv[6] = "tasks.txt";
	cout << "***DEBUG MODE***" << endl;
#endif
	
	if ( argc < 4 ) {
		cout << "Latin square conseq solver: [N] [rows_count, including fixed 1st row] [K] <problem_type> <values_file> <tasks_file>" << endl;
	    return 1;
	}

	cout << "Start solve in latin mode" << endl;

	ls.N          = atoi( argv[1] );
	ls.rows_count = atoi( argv[2] );
	ls.K          = atoi( argv[3] );
	cout << "argc " << argc << endl;
	if ( argc > 4 ) {
		str = argv[4];
		cout << "argv[4] " << str << endl;
		ls.problem_type = str;
	}
	if ( argc > 5 ) {
		IsValuesFile = true;
		vector<int> vec_int;
		cout << "IsValuesFile " << IsValuesFile << endl;
		values_file.open( argv[5] );
		while ( getline( values_file, str ) ) {
			sstream << str;
			while ( sstream >> int_val )
				vec_int.push_back( int_val );
			known_values_vec.push_back( vec_int );
			vec_int.clear();
			sstream.clear(); sstream.str("");
		}
		values_file.close();
	}
	if ( argc > 6 ) { // get number of tasks for creation
		IsTasksFile = true;
		cout << "IsTasksFile " << IsTasksFile << endl;
		tasks_file.open( argv[6] );
		while ( getline( tasks_file, str ) )
			tasks_vec.push_back( atoi( str.c_str()) );
		tasks_file.close();
		// use max number of task for creation
		ls.max_values_len = tasks_vec[ tasks_vec.size()-1 ] * PROBLEMS_IN_WU;
	}

	if ( !IsTasksFile )
		ls.max_values_len   = MAX_VALUES_LEN;
	ls.max_nof_restarts = MAX_NOF_RESTARTS;
	ls.skip_values      = 0;
	ls.verbosity        = 1;
	ls.solver_type      = 4;

	cout << "N "                << ls.N                << endl;
	cout << "rows_count "       << ls.rows_count       << endl;
	cout << "K "                << ls.K                << endl;
	cout << "problem_type "     << ls.problem_type     << endl;
	cout << "IsTasksFile "      << IsTasksFile         << endl;
	cout << "max_nof_restarts " << ls.max_nof_restarts << endl;
	cout << "max_values_len "   << ls.max_values_len   << endl;
	cout << "skip_values "      << ls.skip_values      << endl;
	cout << "verbosity "        << ls.verbosity        << endl;
	cout << "solver_type "      << ls.solver_type      << endl;

	if ( ls.problem_type == "diag" ) {
		ls.cnf_array.resize( sizeof(diag10_2_cnf_array)  / sizeof(diag10_2_cnf_array[0]) );
		for ( unsigned i = 0; i < ls.cnf_array.size(); i++ ) 
			ls.cnf_array[i] = diag10_2_cnf_array[i];
	}

	if ( known_values_vec.size() > 0 )
		ls.known_values_vec = known_values_vec;
	
	ls.MakeLatinValues( );

	if ( IsTasksFile ) { // get PROBLEMS_IN_WU assumptions for every task
		cout << "Start IsTasksFile mode" << endl;
		ofstream ofile( "out_assumptions.txt" );
		int assumption_index;
		for ( unsigned i = 0; i < tasks_vec.size(); i++ )
			for ( unsigned j = 0; j < PROBLEMS_IN_WU; j++ ) {
				assumption_index = (tasks_vec[i]-1) * PROBLEMS_IN_WU + j;
				//ofile << "assumption_index" << assumption_index << endl; 
				ofile << "c 1 12 23 34 45 56 67 78 89 100 "; 
				for ( unsigned j2=0; j2 < ls.positive_literals[assumption_index].size(); j2++ ) {
					ofile << ls.positive_literals[assumption_index][j2];
					if ( j2 < ls.positive_literals[assumption_index].size()-1 )
						ofile << " ";
				}
				ofile << endl;
			}
		ofile.close();
	}
	
	cout << "Make_Permutations_n_k() done" << endl;
	ls.SolveLatinProblems( );
	
	ofstream outfile;
	outfile.open( "outfile", ios_base :: out );
	outfile << ls.all_answers.rdbuf();
	outfile.close();

	return 0;
}