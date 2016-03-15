#include "latin_squares.h"

#include <vector>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "addit_func.h"

#ifdef _DEBUG
	const unsigned MAX_VALUE_COUNT = 100;
#else
	const unsigned MAX_VALUE_COUNT = 5000000;
#endif

using namespace Addit_func;

const int CHECKPOINT_EVERY_SEC = 5;

// --------------------- //
latin_square :: latin_square() :
	N( 0 ),
	K( 0 ),
	rows_count( 0 ),
	diag_elements ( 0 ),
	IsSATFinded( false ),
	verbosity( 0 ),
	problems_solved( 0 ),
	sat_solved( 0 ), 
	unsat_solved( 0 ),
	interrupted( 0 ),
	all_problems( 0 ),
	min_time_sat( 0 ),
	max_time_sat( 0 ), 
	mid_time_sat( 0 ), 
	min_time_unsat( 0 ),
	max_time_unsat( 0 ),
	mid_time_unsat( 0 ),
	min_time_inter( 0 ),
	max_time_inter( 0 ),
	mid_time_inter( 0 ),
	skip_values( 0 ),
	values_checked( 0 ),
	max_nof_restarts( 10 ),
	total_sat_time( 0 ),
	total_unsat_time( 0 ),
	total_inter_time( 0 ),
	max_values_len(MAX_VALUES_LEN),
	out_file_name( "out.txt" ),
	sat_file_name( "sat_sets.txt" ),
	solver_type( 0 ),
	final_values_index( 0 ),
	problem_type ( "diag" ),
	ls_system_rank ( 2 )
{ }

bool latin_square :: ReadLiteralsFromFile( string &input_path, string &error_msg )
{
	string str;
	ifstream ifile( input_path.c_str() );
	getline(ifile, str);
	if ( str.find("diag_start") != string::npos )
		problem_type = "diag";
	else
	{
		error_msg += "impossible head of file " + str;
        return false;
    }
	
	// read and parse string from from file
	string word;
	int val;
	vector<int> cur_vec;
	stringstream sstream;
	while (getline(ifile, str)) {
		sstream << str;
		if ((str.find("p cnf") != std::string::npos) ||
			(str.find("diag_end") != std::string::npos))
		{
			break;
		}
		sstream >> word;
		if (word == "c") { // read literals
			while (sstream >> word) {
				istringstream(word) >> val;
				cur_vec.push_back(val);
			}
			positive_literals.push_back(cur_vec);
			cur_vec.clear();
		}
		sstream.str(""); sstream.clear();
	}

	ifile.close();

	return true;
}

// --------------------- //
void latin_square :: WriteCurrentState( ofstream &out_file, double current_time )
{
// Wtrie current state of latin squares processing
	stringstream sstream;
	//sstream.width( 10 ); sstream.precision( 5 );
	sstream << "*** Current latin state" << endl;
	sstream << "current_time "    << current_time << endl;
	sstream << "positive_literals.size()" << positive_literals.size() << endl;
	sstream << "solved "          << problems_solved << " / " << all_problems << endl;
	sstream << "sat solved "      << sat_solved << endl;
	sstream << "unsat solved "    << unsat_solved << endl;
	sstream << "interrupted "     << interrupted << endl;
	sstream << "min_time_sat "    << min_time_sat << endl;
	sstream << "max_time_sat "    << max_time_sat << endl;
	sstream << "mid_time_sat "    << mid_time_sat << endl;
	sstream << "min_time_unsat "  << min_time_unsat << endl;
	sstream << "max_time_unsat "  << max_time_unsat << endl;
	sstream << "mid_time_unsat "  << mid_time_unsat << endl;
	sstream << "min_time_inter "  << min_time_inter << endl;
	sstream << "max_time_inter "  << max_time_inter << endl;
	sstream << "mid_time_inter "  << mid_time_inter << endl;
	sstream << "***" << endl;

	if ( current_time ) {
		out_file.open( out_file_name.c_str(), ios_base :: out | ios_base :: app );
		out_file << sstream.rdbuf();
		out_file.close();
	}
	cout << sstream.str();
}

// --------------------- //
bool latin_square :: SolveOneProblem( Solver *&S, vector< vector<int> > :: iterator &positive_literals_it, clock_t clk_start )
{
// solve 1 SAT problem encoding latin pair searching problem 
	stringstream sstream;
	vector<int> :: iterator it2;
	unsigned nassigns = 0;
	int st;
	clock_t clk;
	double solving_time = 0;
	vec<Lit> dummy;
	int parsed_lit, var;
	
	for ( it2 = positive_literals_it->begin(); it2 != positive_literals_it->end(); it2++ ) {
		parsed_lit = *it2;
		var = parsed_lit - 1;
		dummy.push( mkLit(var) );
		nassigns++;
	}

	if ( !S->simplify() ) {
        printf("Solved by unit propagation\n");
        printf("\n");
        printf("UNSATISFIABLE\n");
    }

	clk = clock();
    lbool ret = S->solveLimited( dummy );
	solving_time = (double)(clock( ) - clk)/(double)(CLOCKS_PER_SEC);

	if ( ret == l_True )
		st = 1;
	else if ( ret == l_False )
		st = 0;
	else
		st = -1; // l_Undef
	
	problems_solved++;

	if ( st == 1 ) {
		sat_solved++;
		total_sat_time += solving_time;
		mid_time_sat = total_sat_time / (double)sat_solved;
		if ( sat_solved == 1 )
			min_time_sat = max_time_sat = solving_time;
		if ( solving_time < min_time_sat )
			min_time_sat = solving_time;
		if ( solving_time > max_time_sat )
			max_time_sat = solving_time;
		//sat_answers << " SAT ";
		all_answers << " SAT ";
		for ( int k = 0; k < S->model.size() ; k++ ) {
			if ( S->model[k] == l_True ) {
				sstream << k+1 << " ";
				all_answers << "1";
			}
			else  {
				sstream << "-" << k+1 << " ";
				all_answers << "0";
			}
		}
		all_answers << " time " << solving_time << " s" << endl;
		sstream << endl;
		if ( verbosity > 0 ) {
			sat_file.open( sat_file_name.c_str(), ios_base :: out | ios_base :: app );
			sat_file << sstream.rdbuf();
			sat_file.close();
		}
		sstream.clear(); sstream.str("");
	}
	else if ( st == 0 ) {
		unsat_solved++;
		total_unsat_time += solving_time;
		mid_time_unsat = total_unsat_time / (double)unsat_solved;
		all_answers << " UNSAT ";
		all_answers << "time " << solving_time << " s" << endl;
		if ( unsat_solved == 1 )
			min_time_unsat = max_time_unsat = solving_time;
		if ( solving_time < min_time_unsat )
			min_time_unsat = solving_time;
		if ( solving_time > max_time_unsat )
			max_time_unsat = solving_time;
	}
	else if ( st == -1 ) {
		interrupted++;
		interrupted_problems.push_back( *positive_literals_it );
		total_inter_time += solving_time;
		mid_time_inter = total_inter_time / (double)interrupted;
		all_answers << " INTERRUPTED ";
		all_answers << "time " << solving_time << " s" << endl;
		if ( interrupted == 1 )
			min_time_inter = max_time_inter = solving_time;
		if ( solving_time < min_time_inter )
			min_time_inter = solving_time;
		if ( solving_time > max_time_inter )
			max_time_inter = solving_time;
		if ( interrupted % 10000 == 0 )
			cout << "interrupted " << interrupted << endl;
	}

	return true;
}

void latin_square :: SolveLatinProblems( )
{
	cout << "Start FindLatinSquares()" << endl;
	all_problems = (unsigned)positive_literals.size();
	cout << "all_problems " << all_problems << endl; 

	if ( verbosity > 0 ) {
		out_file.open( out_file_name.c_str(), ios_base :: out | ios_base :: app );
		sat_file.open( sat_file_name.c_str(), ios_base :: out );
		sat_file.close();
		out_file << "CNFfile " << CNFfile << " N " << N << " K " << K << " rows_count " << rows_count << endl;
		out_file.close();
	}

	clock_t clk_start = clock( );
	vector< vector<int> > :: iterator positive_literals_it;
	double current_time = 0, time_last_checkpoint = 0;

	minisat22_wrapper m22_wrapper;

	Problem cnf;
	m22_wrapper.parse_DIMACS_from_inc( cnf_array, cnf );
	Solver *S;
	S = new Solver();
	S->addProblem( cnf ); // add initial CNF every time
	S->max_nof_restarts = max_nof_restarts;

	if ( problem_type == "diag" ) // disable reduceDB() each restart
		S->problem_type = problem_type;
	
	for ( positive_literals_it = positive_literals.begin(); positive_literals_it != positive_literals.end(); positive_literals_it++ ) {
		if ( !SolveOneProblem( S, positive_literals_it, clk_start ) ) {
			cout << "Interrupted in client.cpp after SolveOneProblem() returned false" << endl;
			break;
		}
		if ( verbosity > 0 ) {
			current_time = (double)(clock( ) - clk_start)/(double)(CLOCKS_PER_SEC);
			if ( current_time > time_last_checkpoint + CHECKPOINT_EVERY_SEC ) {
				time_last_checkpoint = current_time;
				WriteCurrentState( out_file, current_time );
			}
		}
		//delete S;
	}
	if ( verbosity > 0 ) {
		double current_time = (double)(clock( ) - clk_start)/(double)(CLOCKS_PER_SEC);
		WriteCurrentState( out_file, current_time );
	}
	delete S;
}

bool latin_square :: CompareWithFirstRow( vector<char> vec, unsigned row_index, unsigned columns_count )
{
	// check originality of every value in a row
	for ( unsigned int i = 0; i < vec.size() - 1; i++ ) {
		for ( unsigned int j = i + 1; j < vec.size(); j++ ) {
			if ( vec[i] == vec[j] )
				return false;
		}
	}
	// check columns with 1st known row 1..N
	for ( unsigned i = 0; i < columns_count; i++ ) {
		if ( vec[i] == (char)i + '0' )
			return false;
	}

	if ( problem_type != "diag" )
		return true;
	
	// check if diag cell contains '0' (1st known row has one already)
	if ( ( row_index < vec.size() ) && ( vec[row_index] == '0' ) )
		return false;
	// check if secondary diag cell contains 'N' (1st known row has one already)
	if ( ( N - 1 - row_index < vec.size() ) && ( vec[N - 1 - row_index] == char(N-1) + '0' ) )
		return false;
	
	return true;
}

bool latin_square :: CheckValue( vector<char> cur_vec, int columns_count )
{
// check rows with each other except 1st one, which was checked already
	if ( rows_count == 2 )
		return true;

	int index1, index2;
	for ( int column_index = 0; column_index < columns_count; column_index++ ) {
		for ( int j = 0; j < rows_count - 2; j++ )  {
			for ( int j2 = j + 1; j2 < rows_count-1; j2++ )  {
				index1 = j*columns_count  + column_index;
				index2 = j2*columns_count + column_index;
				if  ( ( (unsigned)index1 < cur_vec.size() ) &&
					  ( (unsigned)index2 < cur_vec.size() ) &&
					  ( cur_vec[index1] == cur_vec[index2] )
					  )
					return false;
			}
		}
	}

	if ( problem_type != "diag" )
		return true;

	// check main and secondary diag. Skip checking 1st row, it was checked already
	for ( int i = 0; i < rows_count - 2; i++ ) {
		for ( int j = i + 1; j < rows_count-1; j++ ) {
			// check main diag. ( N-columns_count  ) is additinal cells of 1st row
			index1 = i*columns_count + i + 1;
			index2 = j*columns_count + j + 1;
			if ( ( i+1    < columns_count  ) && // if index still in current row
				 ( j+1    < columns_count  ) && // if index still in current row
				 ( cur_vec[index1] == cur_vec[index2] )
				 )
				return false;
			// check secondary diag
			index1 = N - i - 2;
			index2 = N - j - 2;
			if ( ( index1 < columns_count ) && // if index still in current row
				 ( index2 < columns_count ) && // if index still in current row
				 ( cur_vec[i*columns_count + index1] == cur_vec[j*columns_count + index2] )
				 )
				return false;
		}
	}

	return true;
}

void latin_square :: FindAdditValues( vector<char> cur_vec, unsigned row_index, vector<char> &add_values )
{
	int count, index;
	for ( int i = 0; i < N; i++ ) {
		count = 0;
		for ( int j = 0; j < K; j++ ) {
			index = (N-K) + K*row_index + j;
			if ( ((unsigned)index < cur_vec.size() ) &&
			     ( cur_vec[index] != char(i) + '0' )
			     )
				count++;
			if ( count == K )
				add_values.push_back( char(i) + '0' );
		}
	}
}

bool latin_square :: IsPossibleValue( vector<char> cur_vec )
{
	if ( !CheckValue( cur_vec, K ) ) // check rows with each other
		return false;
	
	if ( N-K >= 4 ) // too much to check. possible values will exist with high probability
		return true;
	
// each row can be filled up by many variants
// if no variants are possible, then skip current value
	unsigned variants_count_one_row = 1;
	for ( int i=2; i <= N-K; i++ )
		variants_count_one_row *= i;
	
	// find all possible additional values via permutations
	vector< vector<char> > known_rows;
	vector< vector< vector<char> > > add_values;
	known_rows.resize( rows_count - 1 );
	add_values.resize( rows_count - 1 );
	//for ( auto &x : add_values )
	//	x.resize( 1 );
	for( unsigned i=0; i < add_values.size(); ++i )
		add_values[i].resize( 1 );
	// fill first possible value for every row
	for ( unsigned i=0; i < known_rows.size(); i++ ) {
		known_rows[i].resize( K );
		for ( unsigned j=0; j < known_rows[i].size(); j++ )
			known_rows[i][j] = cur_vec[i*K + j];
		for ( int j=0; j < N; j++ ) {
			if ( find(known_rows[i].begin(), known_rows[i].end(), (char)j + '0') == known_rows[i].end() )
				add_values[i][0].push_back( (char)j + '0' );
		}
	}
	// fill all possible values for every row
	vector<char> vec_ch;
	vector< char > :: iterator ch_it;
	for ( unsigned i=0; i < add_values.size(); i++ ) {
		vec_ch = add_values[i][0];
		while ( next_permutation( vec_ch.begin(), vec_ch.end() ) )
			add_values[i].push_back( vec_ch );
	}
	vector<int> index_arr;
	vector< vector<char> > cur_cartesian;
	vector< vector<char> > cur_row_vecs;
	cur_row_vecs.resize( rows_count-1 );
	for( unsigned i=0; i < cur_row_vecs.size(); ++i ){
		cur_row_vecs[i] = known_rows[i];
		cur_row_vecs[i].resize(N);
	}
	unsigned extd_cur_vec_added;
	vector<char> extd_cur_vec;
	bool IsSkip;
	extd_cur_vec.resize( (rows_count-1)*N );
	while( next_cartesian( add_values, index_arr, cur_cartesian ) ) {
		extd_cur_vec_added = 0;
		for ( unsigned i=0; i < known_rows.size(); ++i ) {
			copy( known_rows[i].begin(), known_rows[i].end(), extd_cur_vec.begin() + extd_cur_vec_added );
			extd_cur_vec_added += (unsigned)known_rows[i].size();
			copy( cur_cartesian[i].begin(), cur_cartesian[i].end(), extd_cur_vec.begin() + extd_cur_vec_added );
			extd_cur_vec_added += (unsigned)cur_cartesian[i].size();
		}

		IsSkip = false;
		for ( int i = 0; i < rows_count-1; i++ ) { // check with 1st fixed row 1..N
			for ( int j = N-K+1; j < N; j++ )
				cur_row_vecs[i][j] = extd_cur_vec[N*i + j];
			if ( !CompareWithFirstRow( cur_row_vecs[i], i+1, N ) ) {
				IsSkip = true;
				break;
			}
		}
		if ( IsSkip ) continue;

		if ( CheckValue( extd_cur_vec, N ) ) // check value itself
			return true; // possible variant
	}

	return false;
}

void latin_square :: MakePositiveLiterals()
{
// formula p*n^3 + i*n^2 + j*n + z (p - number of square, i - num of row, j - num of column, z - value )
	stringstream sstream;
	int tmp;
	positive_literals.resize( final_values.size() );
	for ( unsigned i = 0; i < positive_literals.size(); i++ )
		positive_literals[i].resize( final_values[i].size() );
	// rows of 1st square
	for ( unsigned i = 0; i < final_values.size(); i++ )
		for ( int row_index = 1; row_index < rows_count; row_index++ )
			for ( int j = 0; j < K; j++ ) {
				sstream << final_values[i][(row_index-1)*K + j];
				sstream >> tmp;
				positive_literals[i][(row_index-1)*K + j] = row_index*N*N + N*j + (tmp + 1);
				sstream.clear(); sstream.str("");
			}
}

void latin_square :: Show_Values()
{
	/*for ( unsigned i = 0; i < final_values.size(); i++ ) {
		for ( unsigned j = 0; j < final_values[i].size(); j++ )
			cout << final_values[i][j] << " ";
		cout << endl;
		for ( unsigned j = 0; j < final_values[i].size(); j++ )
			cout << positive_literals[i][j] << " ";
		cout << endl;
	}*/
	cout << "Showing first 10 vectors of positive_literals" << endl;
	for ( unsigned int i = 0; i < 10; i++ ) {
		for ( unsigned j = 0; j < positive_literals[i].size(); j++ )
			cout << positive_literals[i][j] << " ";
		cout << endl;
	}
}

// ---------- //
void latin_square :: MakeLatinValues( )
{
	vector< vector< vector<char> > > row_values;
	vector<char> first_row, cur_row;
	vector< vector<int> > permutations;
	
	row_values.resize( rows_count-1 );
	if ( known_values_vec.size() > row_values.size() - 1 )
		known_values_vec.resize( row_values.size() - 1 );
	first_row.resize( N ); // N - order of squares
	for ( int i = 0; i < N; i++ )
		first_row[i] = (char)i + '0';
	cur_row.resize(K);
	
	vector< vector<int> > :: iterator known_values_it, perm_it;
	// Make possible values for every row. 1st known row 1..N has row_index 0 and
	int row_index = 0;
	for ( known_values_it = known_values_vec.begin(); known_values_it != known_values_vec.end(); ++known_values_it ) {
		if ( (*known_values_it).size() < (unsigned)K ) {
			cerr << "Error. (*known_values_it).size() < K" << endl;
			exit(1);
		}
		cout << "known_values_vec" << endl;
		for ( int i=0; i < K; ++i )
			cout << (*known_values_it)[i] << " ";
		cout << endl;
		row_values[row_index].resize(1);
		for ( int i=0; i < K; ++i )
			row_values[row_index][0].push_back( (char)(*known_values_it)[i] + '0' );
		row_index++;
	}
	
	if ( row_index > rows_count-2 ) {
		std::cerr << "row_index > rows_count-2" << ::endl;
		exit(1);
	}
	
	MakePermutations(N, K, permutations); // make pemutations once and use it for every row
	while ( row_index < rows_count-1 ) {
		for ( perm_it=permutations.begin(); perm_it != permutations.end(); ++perm_it ) {
			for ( unsigned j=0; j < (*perm_it).size(); ++j )
				cur_row[j] = (char)(*perm_it)[j] + '0'; // convert to char
			if ( CompareWithFirstRow( cur_row, row_index+1, K ) )
				row_values[row_index].push_back( cur_row );
		}
		cout << "row_values " << row_index << " size " << row_values[row_index].size() << endl;
		cout << "***Permutation() # " << row_index << " done" << endl;
		row_index++;
	}
	
	cout << "short values created" << endl;
	stringstream sstream;
	unsigned long long max_final_values_size = 1;
	for ( int i=0; i<rows_count-1; ++i )
		max_final_values_size *= (unsigned long long)row_values[i].size();
	cout << "max_final_values_size " << max_final_values_size << endl;
	cout << "max_values_len " << max_values_len << endl;
	unsigned long long final_values_size;
	final_values_size = ( max_final_values_size < max_values_len ) ? max_final_values_size : max_values_len;
	cout << "max final values size " << final_values_size << endl;
	final_values.reserve( (unsigned)final_values_size );
	
	vector< vector<char> > row_set;
	vector<char> cur_vec;
	cur_vec.resize( (rows_count-1)*K );
	values_checked = 0;
	unsigned long long impossible_count=0;
	vector<int> index_arr;
	unsigned k;
	while( next_cartesian( row_values, index_arr, row_set ) ) {
		values_checked++;
		if ( ( skip_values ) && ( values_checked <= skip_values ) ) // skip some values
			continue;
		//for( auto &x : row_set )
		//	for( auto &y : x )
		//		cur_vec[i++] = y;
		k=0;
		for( unsigned i=0; i < row_set.size(); ++i )
			for( unsigned j=0; j < row_set[i].size(); ++j )
				cur_vec[k++] = row_set[i][j];
		if ( IsPossibleValue( cur_vec ) ) {
			final_values.push_back( cur_vec );
			if ( final_values.size() == final_values_size ) {
				cout << "final_values.size() == final_values_size. break in next_cartesian()" << endl;
				break;
			}
		}
		else
			impossible_count++;
		if( values_checked % 1000000 == 0 ) {
			cout << "values_checked " << values_checked << endl;
			cout << "impossible_count " << impossible_count << endl;
		}
	}

	row_values.clear();
	cout << "values_checked " << values_checked << endl;
	cout << "impossible_count " << impossible_count << endl;
	cout << "final_values.size() " << final_values.size() << endl;
	
	cout << "before MakePositiveLiterals()" << endl;
	MakePositiveLiterals();
	final_values.clear();
	cout << "positive_literals created" << endl;

	Show_Values();

	sstream << "positive_literals.size() " << positive_literals.size() << endl << endl;
	ofstream out_file;
	out_file.open( out_file_name.c_str(), ios_base :: out );
	out_file << sstream.rdbuf();
	out_file.close();
	cout << sstream.str();
}

void latin_square::makeDiagonalElementsValues()
{
	// make all possible values of diag_elements of first cells of the main and secondary diagonal
	std::cout << "Start makeDiagonalElementsValues()" << std::endl;
	vector< vector<char> > values_main_diag, values_secondary_diag, permutations;
	vector< vector<int> > permutations_int;
	vector<char> cur_value;
	int pemutation_size = 9;
	MakePermutations(N, pemutation_size, permutations_int);
	permutations.resize(permutations_int.size());
	for (unsigned i = 0; i < permutations_int.size(); i++) {
		permutations[i].resize(permutations_int[i].size());
		for (unsigned j = 0; j < permutations_int[i].size(); j++)
			permutations[i][j] = char(permutations_int[i][j]) + '0';
	}
	permutations_int.clear();
	
	std::cout << "permutations.size() " << permutations.size() << std::endl;
	vector<char> fixed_first_row;
	for (int i = 0; i < N; i++)
		fixed_first_row.push_back(char(i) + '0');
	bool isPossibleMainDiag, isPossibleSecondaryDiag, isPossibleValue;
	values_checked = 0;
	values_main_diag.reserve(permutations.size());
	values_secondary_diag.reserve(permutations.size());
	unsigned values_main_diag_count = 0, values_secondary_diag_count = 0;
	unsigned long long values_count = 0;
	
	for (unsigned i = 0; i < permutations.size(); i++) {
		// compare with cell from the fixed first row
		if (find(permutations[i].begin(), permutations[i].end(), fixed_first_row[0]) != permutations[i].end())
			continue;
		// check condition of the main diagonal
		isPossibleMainDiag = true;
		for (unsigned j = 0; j < permutations[i].size(); j++) {
			if (permutations[i][j] == fixed_first_row[j + 1]) {
				isPossibleMainDiag = false;
				break;
			}
		}
		if (isPossibleMainDiag) {
			values_main_diag.push_back(permutations[i]);
			values_main_diag_count++;
		}
	}
	
	values_main_diag.resize(values_main_diag_count);
	std::cout << "values_main_diag.size() " << values_main_diag.size() << std::endl;
	std::stringstream sstream;
	
	if (diag_elements == pemutation_size)
		final_values = values_main_diag;
	else {
		for (unsigned i = 0; i < permutations.size(); i++) {
			// compare with cell from the fixed first row
			if (find(permutations[i].begin(), permutations[i].end(), fixed_first_row[fixed_first_row.size() - 1]) != permutations[i].end())
				continue;
			// check condition of the secondary diagonal
			isPossibleSecondaryDiag = true;
			for (unsigned j = 0; j < permutations[i].size(); j++) {
				if (permutations[i][j] == fixed_first_row[j]) {
					isPossibleSecondaryDiag = false;
					break;
				}
			}
			if (isPossibleSecondaryDiag) {
				values_secondary_diag.push_back(permutations[i]);
				values_secondary_diag_count++;
			}
		}
		
		values_secondary_diag.resize(values_secondary_diag_count);
		std::cout << "values_secondary_diag.size() " << values_secondary_diag.size() << std::endl;
		bool isValueBreak = false;
		
		for (unsigned i = 0; i < values_main_diag.size(); i++) {
			if (isValueBreak)
				break;
			for (unsigned j = 0; j < values_secondary_diag.size(); j++) {
				values_checked++;
				if ((skip_values) && (values_checked <= skip_values)) // skip some values
					continue;
				isPossibleValue = true;
				// compare main and secondary diagonals by N-1 rows (in rows 0 they don't intersect)
				for (int row_index = 1; row_index < N; row_index++) {
					if (values_main_diag[i][row_index-1] == values_secondary_diag[j][N - row_index - 1]) {
						isPossibleValue = false;
						break;
					}
				}
				if (isPossibleValue) {
					// compare main and secondary diagonals by N-2 columns (in column 0 and N-1 they don't intersect)
					for (int col_index = 1; col_index < N-1; col_index++) {
						if (values_main_diag[i][col_index - 1] == values_secondary_diag[j][col_index]) {
							isPossibleValue = false;
							break;
						}
					}
				}
				
				if (isPossibleValue) {
					//values_count++;
					//if (values_count % 1000000 == 0)
					//	std::cout << values_count << std::endl;
					cur_value = values_main_diag[i];
					for (unsigned j2 = 0; j2 < values_secondary_diag[j].size(); j2++ )
						cur_value.push_back(values_secondary_diag[j][j2]);
					final_values.push_back(cur_value);
					//if (final_values.size() % 1000000 == 0)
					//	std::cout << final_values.size() << std::endl;
					if (final_values.size() == max_values_len) {
						cout << "final_values.size() == max_values_len " << max_values_len << ". Break." << endl;
						isValueBreak = true;
						break;
					}
				}
			}
		}
	}
	
	std::cout << "final_values.size() " << final_values.size() << std::endl;
	//std::cout << "values_count " << values_count << std::endl;
	makeDiagonalElementsPositiveLiterals();
}

void latin_square::makeDiagonalElementsPositiveLiterals()
{
	// formula p*n^3 + i*n^2 + j*n + z (p - number of square, i - num of row, j - num of column, z - value )
	positive_literals.resize(final_values.size());
	unsigned row_index, column_index, secondary_diag_index;
	int diagonal_size = N - rows_count;
	int tmp;
	// add literals for the main diagonal
	for (unsigned i = 0; i < final_values.size(); i++)
		for (int diag_element_index = 0; diag_element_index < diagonal_size; diag_element_index++) {
			row_index = column_index = rows_count + diag_element_index;
			tmp = (int)(final_values[i][diag_element_index] - '0');
			positive_literals[i].push_back( row_index*N*N + column_index*N + tmp + 1 );
		}
	if (diag_elements > diagonal_size) {
		// add literals for the secondary diagonal
		for (unsigned i = 0; i < final_values.size(); i++) {
			for (int diag_element_index = diagonal_size; diag_element_index < diagonal_size*2; diag_element_index++) {
				secondary_diag_index = diag_element_index - diagonal_size;
				row_index = N - 1 - secondary_diag_index;
				column_index = secondary_diag_index;
				tmp = (int)(final_values[i][diag_element_index] - '0');
				positive_literals[i].push_back( row_index*N*N + column_index*N + tmp + 1 );
			}
			sort(positive_literals[i].begin(), positive_literals[i].end());
		}
	}
}

void latin_square::makeCnfsFromPositiveLiterals(std::vector<std::string> &base_cnf_vec, int cur_positive_index)
{
	std::string cur_cnf_name, str;
	std::stringstream sstream, base_sstream;
	std::ofstream ofile;
	int tmp;
	
	for (unsigned base_cnf_vec_index = 0; base_cnf_vec_index < base_cnf_vec.size(); base_cnf_vec_index++) {
		std::ifstream ifile(base_cnf_vec[base_cnf_vec_index]);
		while (getline(ifile, str))
			base_sstream << str << std::endl;
		ifile.close(); ifile.clear();

		for (unsigned positive_literals_index = 0;
			positive_literals_index < positive_literals.size();
			positive_literals_index++)
		{
			sstream << cur_positive_index;
			cur_cnf_name = base_cnf_vec[base_cnf_vec_index] + "_" + sstream.str();
			sstream.clear(); sstream.str("");
			cur_cnf_name += ".cnf";
			ofile.open(cur_cnf_name.c_str(), ios::out);
			ofile << base_sstream.rdbuf();

			// fixed normalized rows for every DLS from a system
			for (unsigned j = 0; j < ls_system_rank; j++)
				for (int j2 = 0; j2 < N; j2++) {
					tmp = j*N*N*N + j2*N + j2 + 1; // row index == 0 here
					ofile << tmp << " 0" << std::endl;
				}
			for (unsigned j = 0; j < positive_literals[positive_literals_index].size(); j++) {
				ofile << positive_literals[positive_literals_index][j] << " 0";
				if (j != positive_literals[positive_literals_index].size() - 1)
					ofile << std::endl;
			}
			ofile.close(); ofile.clear();
		}
		
		base_sstream.str("");
		base_sstream.clear();
	}
}

void latin_square::makePositiveLiteralsFromKnownDls( dls known_dls )
{
	// make positive literals for diagonals
	positive_literals.resize(1);
	positive_literals[0].resize(0);
	if (diag_elements > 0) {
		final_values.resize(1);
		// main diagonal
		for (int i = rows_count; i < N; i++)
			final_values[0].push_back(known_dls[i][i]);
		for (int i = N-1; i >= rows_count; i--)
			final_values[0].push_back(known_dls[i][N - 1 - i]);
		makeDiagonalElementsPositiveLiterals();
	}
	// add positive literals for rows
	// for every known row
	std::stringstream sstream;
	int tmp;
	for (int row_index = 1; row_index < rows_count; row_index++) 
		for (int col_index = 0; col_index < N; col_index++) {
			sstream << known_dls[row_index][col_index];
			sstream >> tmp;
			positive_literals[0].push_back( row_index*N*N + N*col_index + (tmp + 1) );
			sstream.clear(); sstream.str("");
		}
}

void latin_square::readOdlsPairs(std::string known_podls_file_name)
{
	std::ifstream odls_pairs_file(known_podls_file_name.c_str());
	if (!odls_pairs_file.is_open()) {
		std::cerr << known_podls_file_name << " is not open" << std::endl;
		exit(1);
	}
	std::string str, nonspace_str;
	odls_pair cur_odls_pair;
	while (getline(odls_pairs_file, str)) {
		if ((str.size() <= 1) && (cur_odls_pair.dls_1.size() > 0)) { // if separate string
			odls_pair_vec.push_back(cur_odls_pair); // add current odls pair
			cur_odls_pair.dls_1.clear();
			cur_odls_pair.dls_2.clear();
		}
		else {
			// read two rows of 2 DLS for current pair
			nonspace_str = "";
			for (auto &x : str)
				if ((x != ' ') && (x != '\r'))
					nonspace_str += x;
			if (nonspace_str.size() != 2 * LS_ORDER) {
				std::cerr << "nonspace_str.size() != 2*LS_ORDER" << std::endl;
				std::cerr << nonspace_str.size() << " != " << 2 * LS_ORDER << std::endl;
				exit(1);
			}
			cur_odls_pair.dls_1.push_back(nonspace_str.substr(0, LS_ORDER));
			cur_odls_pair.dls_2.push_back(nonspace_str.substr(LS_ORDER, LS_ORDER));
		}
	}
	odls_pairs_file.close();
	// if there is no separate string at the end of file, add last pair manually
	if (cur_odls_pair.dls_1.size() != 0)
		odls_pair_vec.push_back(cur_odls_pair); // add current odls pair
	
	std::cout << "odls_pair_vec.size() " << odls_pair_vec.size() << std::endl;
	
	std::set<std::string> greece_latin_square;
	std::string cell_plus_cell;
	// check corectness for every pair
	for (auto &x : odls_pair_vec) {
		for (unsigned j1 = 0; j1 < x.dls_1.size(); j1++)
			for (unsigned j2 = 0; j2 < x.dls_1[j1].size(); j2++) {
				cell_plus_cell = x.dls_1[j1][j2];
				cell_plus_cell += x.dls_2[j1][j2];
				greece_latin_square.insert(cell_plus_cell);
			}
		//std::cout << "greece_latin_square.size() " << greece_latin_square.size() << std::endl;
		if (greece_latin_square.size() != x.dls_1.size() * x.dls_1.size()) {
			std::cerr << "greece_latin_square.size() != x.dls_1.size() * x.dls_1.size() " << std::endl;
			std::cerr << greece_latin_square.size() << " != " << x.dls_1.size() * x.dls_1.size() << std::endl;
		}
	}
}

void latin_square::makePseudotriple(odls_pair &orthogonal_pair, dls &new_dls, odls_pseudotriple &pseudotriple)
{
	unsigned cur_first_pair_orthogonal_cells, cur_second_pair_orthogonal_cells;
	std::set<std::string> greece_latin_square1, greece_latin_square2; // for counting orthogonal cells between 2 DLS
	std::string cell_plus_cell;

	for (unsigned j1 = 0; j1< new_dls.size(); j1++)
		for (unsigned j2 = 0; j2 < new_dls[j1].size(); j2++) {
			cell_plus_cell = orthogonal_pair.dls_1[j1][j2];
			cell_plus_cell += new_dls[j1][j2];
			greece_latin_square1.insert(cell_plus_cell);
		}
	cur_first_pair_orthogonal_cells = greece_latin_square1.size();
	/*std::cout << "greece_latin_square1 " << std::endl;
	for ( auto &y : greece_latin_square1 )
	std::cout << y << " ";
	std::cout << std::endl;*/
	for (unsigned j1 = 0; j1 < new_dls.size(); j1++)
		for (unsigned j2 = 0; j2 < new_dls[j1].size(); j2++) {
			cell_plus_cell = orthogonal_pair.dls_2[j1][j2];
			cell_plus_cell += new_dls[j1][j2];
			greece_latin_square2.insert(cell_plus_cell);
		}
	cur_second_pair_orthogonal_cells = greece_latin_square2.size();
	pseudotriple.dls_1 = orthogonal_pair.dls_1;
	pseudotriple.dls_2 = orthogonal_pair.dls_2;
	pseudotriple.dls_3 = new_dls;
	pseudotriple.unique_orthogonal_cells.clear();
	std::set_intersection(greece_latin_square1.begin(), greece_latin_square1.end(),
		greece_latin_square2.begin(), greece_latin_square2.end(),
		std::inserter(pseudotriple.unique_orthogonal_cells, pseudotriple.unique_orthogonal_cells.begin()));
}

void latin_square::makeCnfsFromDls()
{
	readOdlsPairs("ODLS_10_pairs.txt");
	
	std::vector<dls> dls_vec;
	for (unsigned i = 0; i < odls_pair_vec.size(); i++) {
		if ( find(dls_vec.begin(), dls_vec.end(), odls_pair_vec[i].dls_1) == dls_vec.end() )
			dls_vec.push_back(odls_pair_vec[i].dls_1);
		if (find(dls_vec.begin(), dls_vec.end(), odls_pair_vec[i].dls_2) == dls_vec.end())
			dls_vec.push_back(odls_pair_vec[i].dls_2);
	}

	std::cout << "dls_vec.size() " << dls_vec.size() << std::endl;

	if ( max_values_len > dls_vec.size())
		max_values_len = dls_vec.size();

	std::vector<std::string> base_cnf_vec;
	base_cnf_vec.push_back("../../../../Tests/ferrumsat/Latin Square encodings/DLS_10_2_encodings/LSD10_2_pw_ext_bm_ext.cnf");
	base_cnf_vec.push_back("../../../../Tests/ferrumsat/Latin Square encodings/DLS_10_2_encodings/LSD10_2_pw_ext_bn_ext.cnf");
	base_cnf_vec.push_back("../../../../Tests/ferrumsat/Latin Square encodings/DLS_10_2_encodings/LSD10_2_pw_ext_cm_ext.cnf");
	base_cnf_vec.push_back("../../../../Tests/ferrumsat/Latin Square encodings/DLS_10_2_encodings/LSD10_2_pw_ext_pr_ext.cnf");
	base_cnf_vec.push_back("../../../../Tests/ferrumsat/Latin Square encodings/DLS_10_2_encodings/LSD10_2_pw_ext_pw_ext.cnf");
	base_cnf_vec.push_back("../../../../Tests/ferrumsat/Latin Square encodings/DLS_10_2_encodings/LSD10_2_pw_ext_sq_ext.cnf");
	base_cnf_vec.push_back("../../../../Tests/ferrumsat/Latin Square encodings/DLS_10_2_encodings/LSD10_2_pw_naive_pw_naive.cnf");
	
	for (unsigned i = 0; i < max_values_len; i++) {
		makePositiveLiteralsFromKnownDls(dls_vec[i]);
		makeCnfsFromPositiveLiterals(base_cnf_vec, i);
	}
}

void latin_square::makeHtmlData()
{
	std::ifstream ifile("pseudotriple.txt");
	std::string str, tmp_str, cur_dls_row;
	dls cur_dls;
	std::vector<dls> dls_vec;
	std::stringstream sstream;
	while (getline(ifile, str)) {
		if (str.size() < 18) continue;
		sstream << str;
		while (sstream >> tmp_str) {
			cur_dls_row += tmp_str;
			if (cur_dls_row.size() == 10) {
				cur_dls.push_back(cur_dls_row);
				cur_dls_row = "";
			}
			if (cur_dls.size() == 10) {
				dls_vec.push_back(cur_dls);
				cur_dls.clear();
			}
		}
		sstream.clear(); sstream.str("");
	}
	ifile.close();

	std::ofstream html_data_file("html_data.txt");
	for (unsigned i = 1; i < dls_vec.size(); i++) {
		sstream << "<tr>" << std::endl;
		sstream << "<td> " << i << " </td>" << std::endl;
		sstream << "<td>" << std::endl;
		sstream << "<FONT SIZE = -2>" << std::endl;
		for (unsigned j = 0; j < 10; j++) {
			for (unsigned j2 = 0; j2 < 10; j2++)
				sstream << dls_vec[0][j][j2] << " ";
			sstream << "&nbsp;&nbsp ";
			for (unsigned j2 = 0; j2 < 10; j2++)
				sstream << dls_vec[i][j][j2] << " ";
			sstream << "<br>" << std::endl;
		}
		sstream << "</FONT>" << std::endl;
		sstream << "</td>" << std::endl;
		sstream << "</tr>" << std::endl;
	}
	html_data_file << sstream.str();
	html_data_file.close();
}


// normalize LS by 1st row
void normalizeLS(dls &cur_DLS)
{
	char ch;
	dls new_DLS;
	unsigned LS_order = cur_DLS.size();
	new_DLS.resize(cur_DLS.size());
	for (unsigned i = 0; i < LS_order; i++)
		new_DLS[i].resize(LS_order);
	for (unsigned i = 0; i < LS_order; i++) {
		ch = cur_DLS[0][i];
		for (unsigned j = 0; j < LS_order; j++)
			for (unsigned j2 = 0; j2 < LS_order; j2++) {
				if (cur_DLS[j][j2] == ch)
					new_DLS[j][j2] = '0' + i;
			}
	}
	cur_DLS = new_DLS;
}

// Only literals which correspond to 1 are to be made
std::vector<int> makeLiterals(dls cur_dls)
{
	std::vector<int> result_vec;
	unsigned row_index, column_index;
	for (row_index = 1; row_index < cur_dls.size() - 1; row_index++) { // skip row # 0 - it if fixed already
		for (column_index = 0; column_index < cur_dls[0].size() - 1; column_index++)
			result_vec.push_back(1000 * 1 + 100 * row_index + 10 * column_index + (cur_dls[row_index][column_index] - 48) + 1);
	}
	return result_vec;
}

std::vector<dls> getSetUniqueDLS(std::vector<odls_pair> odls_pair_vec)
{
	std::vector<dls> unique_dls_vec;

	bool isInsertRequired;

	for (auto &x : odls_pair_vec) {
		if (x.dls_1[0] != "0123456789") {
			normalizeLS(x.dls_1);
			for (auto &y : x.dls_1) {
				for (auto &y2 : y)
					std::cout << y2 << " ";
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}
		if (x.dls_2[0] != "0123456789") {
			normalizeLS(x.dls_2);
			for (auto &y : x.dls_2) {
				for (auto &y2 : y)
					std::cout << y2 << " ";
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}
		isInsertRequired = true;
		for (auto &y : unique_dls_vec)
			if (y == x.dls_1) {
				isInsertRequired = false;
				break;
			}
		if (isInsertRequired)
			unique_dls_vec.push_back(x.dls_1);
		//
		isInsertRequired = true;
		for (auto &y : unique_dls_vec)
			if (y == x.dls_2) {
				isInsertRequired = false;
				break;
			}
		if (isInsertRequired)
			unique_dls_vec.push_back(x.dls_2);
	}

	return unique_dls_vec;
}