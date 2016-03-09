#ifndef latin_squares_h
#define latin_squares_h

#include <cstdio>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include <set>

#include "minisat22_wrapper.h"
#include "addit_func.h"

typedef std::vector<std::string> dls;

const unsigned LS_ORDER = 10;
const unsigned psuedotriple_char_arr_len = 3 * LS_ORDER * LS_ORDER;

struct odls_pair
{
	dls dls_1;
	dls dls_2;
};

struct odls_pseudotriple
{
	dls dls_1;
	dls dls_2;
	dls dls_3;
	std::set<std::string> unique_orthogonal_cells;
};

const long long MAX_VALUES_LEN = 10000000;

using namespace std;

class latin_square
{
public:
	latin_square();

	int N; // order of squares
	unsigned K; // default number of columns
	unsigned rows_count;
	int diag_elements;
	string problem_type;
	bool IsSATFinded;
	unsigned verbosity;
	int solver_type;
	unsigned long long skip_values;
	unsigned long long values_checked;
	unsigned max_nof_restarts;
	unsigned long long max_values_len;
	vector<int> cnf_array;
	int cnf_array_len;
	vector< vector<int> > known_values_vec;
	vector< vector<int> > positive_literals;
	stringstream all_answers;
	unsigned ls_system_rank;

	bool ReadLiteralsFromFile( string &input_path, string &error_msg );
	void WriteCurrentState( ofstream &out_file, double current_time );
	void Show_Values();

	bool IsPossibleValue( vector<char> cur_vec );
	bool CheckValue( vector<char> cur_vec, unsigned columns_count );
	bool CompareWithFirstRow( vector<char> vec, unsigned row_index, unsigned columns_count );
	void FindAdditValues( vector<char> cur_vec, unsigned row_index, vector<char> &add_values );
	
	void SolveLatinProblems( );
	bool SolveOneProblem( Solver *&S, vector< vector<int> > :: iterator &positive_literals_it, clock_t clk_start );
	
	void MakeLatinValues();
	void MakePositiveLiterals();
	void makeDiagonalElementsValues();
	void makeDiagonalElementsPositiveLiterals();
	void makeCnfsFromPositiveLiterals();
	void makePositiveLiteralsFromKnownDls(dls known_dls);
	void makeCnfsFromDls();
	
	// work with pseudotriples
	void readOdlsPairs(std::string known_podls_file_name);
	void makePseudotriple(odls_pair &orthogonal_pair, dls &new_dls, odls_pseudotriple &pseudotriple);
private:
	std::vector<odls_pair> odls_pair_vec;
	vector< vector<char> > final_values;
	vector< vector<int> > interrupted_problems;
	ofstream out_file;
	ofstream sat_file;
	string CNFfile;
	string out_file_name;
	string sat_file_name;
	unsigned problems_solved;
	unsigned sat_solved;
	unsigned unsat_solved;
	unsigned interrupted;
	unsigned all_problems;
	double total_sat_time;
	double total_unsat_time;
	double total_inter_time;
	double min_time_sat;
	double max_time_sat;
	double mid_time_sat;
	double min_time_unsat;
	double max_time_unsat;
	double mid_time_unsat;
	double min_time_inter;
	double max_time_inter;
	double mid_time_inter;
	string cnf_head_str;
	unsigned final_values_index;
	unsigned out_cnf_global_index;
};

#endif