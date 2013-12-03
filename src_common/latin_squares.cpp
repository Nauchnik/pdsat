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
	max_values_len( 1000000 ),
	out_file_name( "out.txt" ),
	sat_file_name( "sat_sets.txt" ),
	solver_type( 0 ),
	final_values_index( 0 ),
	problem_type ( "inc72" )
{ }

bool latin_square :: ReadLiteralsFromFile( FILE *infile, string &error_msg )
{
	char word1[100], word2[100];
	vector<int> cur_vec;
	int string_count = 0;
	string str1, str2;
	stringstream sstream, total_sstream;

	// read and check head of file
	fscanf( infile, "%s %s", &word1, &word2 );
	str1 = word1;
	str2 = word2;
	if ( !( strcmp( word1, "c" ) ) && !( strcmp( word2, "diag_start" ) ) )
		problem_type = "diag";
	else if ( !( strcmp( word1, "c" ) ) && !( strcmp( word2, "ls_10_3_inc0" ) ) )
		problem_type = "inc0";
	else if ( !( strcmp( word1, "c" ) ) && !( strcmp( word2, "ls_10_3_inc60" ) ) )
		problem_type = "inc60";                               
	else if ( !( strcmp( word1, "c" ) ) && !( strcmp( word2, "ls_10_3_inc70" ) ) )
		problem_type = "inc70";
	else if ( !( strcmp( word1, "c" ) ) && !( strcmp( word2, "ls_10_3_inc70m" ) ) )
		problem_type = "inc70";
	else if ( !( strcmp( word1, "c" ) ) && !( strcmp( word2, "ls_10_3_inc80" ) ) )
		problem_type = "inc80";
	else if ( !( strcmp( word1, "c" ) ) && !( strcmp( word2, "ls_10_3_inc90" ) ) )
		problem_type = "inc90";
	else
	{
		error_msg += "impossible head of file " + str1 + " " + str2;
        return false;
    }
	
	// read from FILE (except head string) to strinstream
	char c = 0;
	while ( c != EOF ) {
		c = fgetc( infile );
		total_sstream << c;
	}
	fclose( infile );

	string word;
	int val;
	while ( !total_sstream.eof() ) {
		total_sstream >> word;
		if ( ( word.find( "cnf" )      != std::string::npos ) ||
		     ( word.find( "diag_end" ) != std::string::npos ) )
		{
			break;
		}
		if ( ( word == "c" ) || ( word == "p" ) ) {
			if ( cur_vec.size() == 0 )
				continue;
			positive_literals.push_back( cur_vec );
			cur_vec.clear();
		}
		else {
			istringstream( word ) >> val;
			cur_vec.push_back( val );
		}
	}

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
		dummy.push( (parsed_lit > 0) ? mkLit(var) : ~mkLit(var) );
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
	all_problems = positive_literals.size();
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

	if ( problem_type == "incall" ) // disable reduceDB() each restart
		S->core_len = S->nVars();
	
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

bool latin_square :: CheckValue( vector<char> cur_vec, unsigned columns_count )
{
// check rows with each other except 1st one, which was checked already
	if ( rows_count == 2 )
		return true;

	unsigned index1, index2;
	for ( unsigned column_index = 0; column_index < columns_count; column_index++ ) {
		for ( unsigned j = 0; j < rows_count - 2; j++ )  {
			for ( unsigned j2 = j + 1; j2 < rows_count-1; j2++ )  {
				index1 = j*columns_count  + column_index;
				index2 = j2*columns_count + column_index;
				if  ( ( index1 < cur_vec.size() ) &&
					  ( index2 < cur_vec.size() ) &&
					  ( cur_vec[index1] == cur_vec[index2] )
					  )
					return false;
			}
		}
	}

	if ( problem_type != "diag" )
		return true;

	// check main and secondary diag. Skip checking 1st row, it was checked already
	for ( unsigned i = 0; i < rows_count - 2; i++ ) {
		for ( unsigned j = i + 1; j < rows_count-1; j++ ) {
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
	unsigned count, index;
	for ( int i = 0; i < N; i++ ) {
		count = 0;
		for ( int j = 0; j < K; j++ ) {
			index = (N-K) + K*row_index + j;
			if ( ( index < cur_vec.size() ) && 
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
	for ( unsigned i=2; i <= N-K; i++ )
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
			extd_cur_vec_added += known_rows[i].size();
			copy( cur_cartesian[i].begin(), cur_cartesian[i].end(), extd_cur_vec.begin() + extd_cur_vec_added );
			extd_cur_vec_added += cur_cartesian[i].size();
		}

		IsSkip = false;
		for ( unsigned i = 0; i < rows_count-1; i++ ) { // check with 1st fixed row 1..N
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

void latin_square :: MakePositiveLiterals( )
{
// formula p*n^3 + i*n^2 + j*n + z (p - number of square, i - num of row, j - num of column, z - value )
	stringstream sstream;
	int tmp;
	positive_literals.resize( final_values.size() );
	for ( unsigned i = 0; i < positive_literals.size(); i++ )
		positive_literals[i].resize( final_values[i].size() );
	// rows of 1st square
	for ( unsigned i = 0; i < final_values.size(); i++ )
		for ( unsigned row_index = 1; row_index < rows_count; row_index++ )
			for ( unsigned j = 0; j < K; j++ ) {
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
	first_row.resize( N ); // N - order of squares
	for ( int i = 0; i < N; i++ )
		first_row[i] = (char)i + '0';
	cur_row.resize(K);
	
	vector< vector<int> > :: iterator known_values_it, perm_it;
	// Make possible values for every row. 1st known row 1..N has row_index 0 and
	unsigned row_index = 0;
	for ( known_values_it = known_values_vec.begin(); known_values_it != known_values_vec.end(); ++known_values_it ) {
		if ( (*known_values_it).size() < K ) {
			cerr << "Error. (*known_values_it).size() < K" << endl;
			exit(1);
		}
		cout << "known_values_vec" << endl;
		for ( unsigned i=0; i < K; ++i )
			cout << (*known_values_it)[i] << " ";
		cout << endl;
		row_values[row_index].resize(1);
		for ( unsigned i=0; i < K; ++i )
			row_values[row_index][0].push_back( (char)(*known_values_it)[i] + '0' );
		row_index++;
	}
	
	if ( row_index > rows_count-2 ) {
		cerr << "row_index > rows_count-2" << endl;
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
	for ( unsigned i=0; i<rows_count-1; ++i )
		max_final_values_size *= (unsigned long long)row_values[i].size();
	cout << "max_final_values_size " << max_final_values_size << endl;
	cout << "max_values_len " << max_values_len << endl;
	unsigned long long final_values_size;
	final_values_size = ( max_final_values_size < max_values_len ) ? max_final_values_size : max_values_len;
	cout << "max final values size " << final_values_size << endl;
	final_values.reserve( final_values_size );
	
	vector< vector<char> > row_set;
	vector<char> cur_vec;
	cur_vec.resize( (rows_count-1)*K );
	values_checked = 0;
	unsigned impossible_count=0;
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