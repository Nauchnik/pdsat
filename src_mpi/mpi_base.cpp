#include "mpi_base.h"

#pragma warning( disable : 4996 )

const int    MAX_WORD_LENGTH			   = 64;
const int    MAX_LINE_LENGTH               = 524288;
const int    MAX_LINE_LENGTH_2             = 8192;
const int    MEDIUM_STRING_LEN             = 4096;
const double TRANSP_COAST                  = 0.000001;
const int    NUM_KEY_BITS                  = 64;

boost::random::mt19937 gen;

//=================================================================================================
// Constructor/Destructor:

MPI_Base :: MPI_Base( ) :
	sat_count            ( 0 ),
	corecount			 ( 0 ),
	solver_type          ( 4 ),
	core_len             ( 0 ),    
	koef_val             ( 64 ),
	schema_type          ( "" ),
	var_count            ( 0 ),
	clause_count         ( 0 ),
	full_mask_var_count  ( 0 ),
	start_activity		 ( 0 ),
	part_mask_var_count  ( 0 ),
	all_tasks_count		 ( 0 ),
	constr_clauses_count ( -1 ),
	obj_clauses_count    ( -1 ),
	obj_vars_count       ( -1 ),
	IsConseq             ( false ),
	IsPB                 ( false ),
	PB_mode				 ( 1 ),
	best_lower_bound	 ( -1 ),
	upper_bound          ( -1 ),
	verbosity			 ( 0 ),
	check_every_conflict ( 2000 ),
	known_point_file_name ( "known_point" ),
	known_assumptions_file_name ( "known_assumptions_0" ),
	base_known_assumptions_file_name ( "known_assumptions" ),
	IsSolveAll           ( false ),
	IsPredict            ( false ),
	max_solving_time     ( 0 ),
	max_nof_restarts     ( 0 ),
	keybit_count         ( 4 ),
	rslos_table_name     ( "" ),
	assumptions_string_count ( 0 ),
	activity_vec_len	 ( 0 )
{
	full_mask  = new unsigned[FULL_MASK_LEN];
	part_mask  = new unsigned[FULL_MASK_LEN];
	mask_value = new unsigned[FULL_MASK_LEN];
	for ( unsigned i = 0; i < FULL_MASK_LEN; i++ )
		full_mask[i] = part_mask[i] = mask_value[i] = 0;

	gen.seed( static_cast<unsigned>(std::time(0)) );
}

MPI_Base :: ~MPI_Base( )
{
	delete[] full_mask;
	delete[] part_mask;
	delete[] mask_value;
}

// Make full_mask and part_mask for sending from order that is set by var choose array
//---------------------------------------------------------
bool MPI_Base :: GetMainMasksFromVarChoose( vector<int> &var_choose_order )
{
	if ( verbosity > 0 ) {
		cout << endl << "GetMainMasksFromVarChoose() start with var_choose_order " << endl;
		for ( unsigned i = 0; i < var_choose_order.size(); i++ )
			cout << var_choose_order[i] << " ";
		cout << endl;
		cout << "part_mask_var_count " << part_mask_var_count << endl;
		cout << endl;
	}

	// init masks every launch
	for ( unsigned i = 0; i < FULL_MASK_LEN; ++i )
		part_mask[i] = full_mask[i] = 0;
	
	unsigned cur_uint_ind, var_index;
	for ( unsigned i = 0; i < var_choose_order.size(); ++i ) {
		var_index = var_choose_order[i] - 1;
		cur_uint_ind = var_index / UINT_LEN;
		full_mask[cur_uint_ind + 1] += 1 << ( var_index % UINT_LEN );
	}

	// get first full_mask_var_count vars from  array var_choose_order
	for ( unsigned i = 0; i < part_mask_var_count; ++i ) {
		var_index = var_choose_order[i] - 1;
		cur_uint_ind = var_index / UINT_LEN;
		part_mask[cur_uint_ind + 1] += 1 << ( var_index % UINT_LEN );
	}

	// full_mask[0] is mask of existing all 32 uint values
	for ( unsigned i = 1; i < FULL_MASK_LEN; ++i ) // fill full_mask[0]
		if ( full_mask[i] ) full_mask[0] += 1 << ( i-1 );
	
	for ( unsigned i = 1; i < FULL_MASK_LEN; ++i ) // fill part_mask[0]
		if ( part_mask[i] ) part_mask[0] += 1 << ( i-1 );

	if ( ( !full_mask[0] ) && ( full_mask_var_count != 0 ) ) {
		cout << "Error. full_mask[0] == 0. full_mask_var_count != 0";
		return false;
	}
	
	if ( verbosity > 1 ) {
		cout << "made part_mask" << endl;
		for( unsigned i = 0; i < FULL_MASK_LEN; ++i )
			cout << part_mask[i] << " ";
		cout << endl;
		unsigned bit_count = 0;
		for ( unsigned j=1; j < FULL_MASK_LEN; ++j )
			bit_count += BitCount( part_mask[j] );
		cout << "part_mask bit_count " << bit_count << endl;
	}
	
	if ( verbosity > 1 )
		cout << "GetMainMasksFromVarChoose() end" << endl;
	
	return true;
}

bool MPI_Base :: MakeAssignsFromFile( int current_task_index, vec< vec<Lit> > &dummy_vec )
{
	if ( verbosity > 0 )
		cout << "MakeAssignsFromFile()" << endl;
	
	ifstream in;
	string str;
	in.open( known_assumptions_file_name.c_str() );
	if ( !in.is_open() ) {
		cerr << "Error. !in.is_open(). file name " << known_assumptions_file_name << endl;
		return false;
	}
	
	if ( var_choose_order.size() == 0 ) {
		cerr << "Error. var_choose_order.size() == 0 " << endl;
		return false;
	}
	
	int cur_var_ind, intval, k = 0;
	// int rslos_num = 1; 
	int basic_batch_size = (int)floor( (double)assumptions_string_count / (double)all_tasks_count );
	// calculate count of bathes with additional size (+1)
	int batch_addit_size_count = assumptions_string_count - basic_batch_size*all_tasks_count;
	int cur_batch_size = basic_batch_size;
	if ( current_task_index < batch_addit_size_count )
		cur_batch_size++;
	// skip unuseful strings
	int previous_tasks_count = current_task_index*basic_batch_size;
	if ( current_task_index < batch_addit_size_count )
		previous_tasks_count += current_task_index; // add some 1 to sum
	else
		previous_tasks_count += batch_addit_size_count;
	
	int strings_passed = 0;
	while ( strings_passed < previous_tasks_count ) {
		getline( in, str );
		strings_passed++;
	}
	if ( verbosity > 0 ) {
		cout << "current_task_index "       << current_task_index << endl;
		cout << "all_tasks_count "          << all_tasks_count << endl;
		cout << "previous_tasks_count "     << previous_tasks_count << endl;
		cout << "assumptions_string_count " << assumptions_string_count << endl;
		cout << "basic_batch_size "         << basic_batch_size  << endl;
		cout << "cur_batch_size "           << cur_batch_size << endl;
		cout << "strings_passed "           << strings_passed  << endl;
	}
	
	dummy_vec.resize( cur_batch_size );
	// reading values from file
	for ( int i=0; i < cur_batch_size; ++i ) {
		if ( !getline( in, str ) ) {
			cerr << "Error. !getline( in, str )" << endl;
			return false;
		}
		if ( str.size() < var_choose_order.size() ) {
			cerr << "Error. str.size() < var_choose_order.size()" << endl;
			return false;
		}
		for ( unsigned j=0; j < var_choose_order.size(); ++j ) {
			cur_var_ind = var_choose_order[j] - 1;
			if ( str[j] == '1' )
				dummy_vec[i].push( mkLit( cur_var_ind ) );
			else 
				dummy_vec[i].push( ~mkLit( cur_var_ind ) );
		}
	}
	in.close();
	return true;
}

bool MPI_Base :: MakeAssignsFromMasks( unsigned *full_mask, 
									   unsigned *part_mask, 
									   unsigned *mask_value,
									   vec< vec<Lit> > &dummy_vec )
{
// for predict with minisat2.2. convert masks to vector of Literals
	unsigned full_mask_var_count = 0, variate_vars_count = 0;
	for ( unsigned i = 1; i < FULL_MASK_LEN; i++ ) {
		variate_vars_count  += BitCount( full_mask[i] ^ part_mask[i] );
		full_mask_var_count += BitCount( full_mask[i] );
	}
	
	// determine count of assumptions and lenght of every one
	int problems_count = 1 << variate_vars_count;
	dummy_vec.resize( problems_count );
	for( int i=0; i < dummy_vec.size(); ++i )
		dummy_vec[i].resize( full_mask_var_count );
	
	unsigned mask, range_mask_ind;
	int range_mask;
	unsigned index;
	int cur_var_ind;
	bool IsPositiveLiteral, IsAddingLiteral;

	for ( int lint = 0; lint < problems_count; lint++ ) {
		range_mask = 1;
		range_mask_ind = 1;
		index = 0;
		for ( unsigned i = 1; i < FULL_MASK_LEN; i++ ) {
			for ( unsigned j = 0; j < UINT_LEN; j++ ) {
				mask = ( 1 << j );
				cur_var_ind = ( i-1 ) * UINT_LEN + j;
				IsAddingLiteral = false;
				if ( part_mask[i] & mask ) { // one common vector send by control process
					IsPositiveLiteral = ( mask_value[i] & mask ) ? true : false;
					IsAddingLiteral = true;
				}
				else if ( full_mask[i] & mask ) { // set of vectors
					IsPositiveLiteral = ( lint & range_mask ) ? true : false;
					range_mask = 1 << range_mask_ind;
					range_mask_ind++;
					IsAddingLiteral = true;
				}
				if ( !IsAddingLiteral ) continue;
				dummy_vec[lint][index++] = IsPositiveLiteral ? mkLit( cur_var_ind ) : ~mkLit( cur_var_ind );
			}
		}
	}
	return true;
}

// Get values for sending using order by var choose array
//---------------------------------------------------------
bool MPI_Base :: GetValuesFromVarChoose( unsigned &part_var_power )
{
	unsigned cur_uint_ind = 0;
	unsigned value_index;
	unsigned mask;
	unsigned mask2;
	for( unsigned lint = 0; lint < part_var_power; lint++ ) {
		value_index = 0;
		for( unsigned i = 1; i < FULL_MASK_LEN; i++ ) {
			for( unsigned j = 0; j < UINT_LEN; j++ ) {
				mask  = ( 1 << j );
				mask2 = 1 << value_index;
				if ( part_mask[i] & mask ) { // if part_mask bit is 1	
					if ( lint & mask2 ) // if bit is 1
						values_arr[lint][i] += mask;
					value_index++;
				}
			} // for( j = 0; j < UINT_LEN; j++ )
		}
	} // for( lint = 0; lint < part_var_power; lint++ )

	// set values_arr[lint][0] and values_arr[lint][i]
	for( unsigned lint = 0; lint < part_var_power; lint++ ) {
		for( unsigned i = 1; i < FULL_MASK_LEN; i++ ) { // fill part_mask[0]
			if ( values_arr[lint][i] )
				values_arr[lint][0] += 1 << ( i-1 );
		}
	}

	return true;
}

//---------------------------------------------------------
bool MPI_Base :: MakeStandardMasks( unsigned &part_var_power )
{		
	if ( !GetMainMasksFromVarChoose( var_choose_order ) ) { 
		cout << "Error in GetMainMasksFromVarChoose" << endl; 
		return false;
	}
	cout << "Correct end of GetMasksFromVarChoose" << endl;
		
	if ( !GetValuesFromVarChoose( part_var_power ) ) { 
		cout << "Error in GetValuesFromVarChoose" << endl; 
		return false; 
	}
	cout << "Correct end of GetValuesFromVarChoose" << endl;
	
	return true;
}

//---------------------------------------------------------
bool MPI_Base :: MakeVarChoose( )
{
// Make array var_choose_order with vars sorted by given rule
	unsigned *var_literal_count_weights;
	unsigned *var_implicant_count_weights;
	double *var_jeroslaw_count_weights;
	
	string str;
	stringstream sstream;
	int val;
	
	// if file with decomposition set exists
	ifstream known_point_file( known_point_file_name.c_str() );
	if ( known_point_file.is_open() ) {
		getline( known_point_file, str );
		sstream << str;
		var_choose_order.resize( 0 );
		while ( sstream >> val )
			var_choose_order.push_back( val );
		full_mask_var_count = var_choose_order.size();
		sstream.str( "" ); sstream.clear( );
		known_point_file.close();
		schema_type = "known_point";
	}
	
	cout << "var_choose_order.size() " << var_choose_order.size() << endl;
	cout << "full_mask_var_count " << full_mask_var_count << endl;

	if ( ( var_choose_order.size() == 0 ) && ( schema_type == "" ) )
		schema_type = "0";
	// if got from known point file or from "c var_set..." then set was made already
	var_choose_order.resize( full_mask_var_count );
	
	if ( schema_type != "" ) {
		int k = 0;
		if ( schema_type == "0" ) { // == bivium_Beginning1
			for ( unsigned i = 0; i < full_mask_var_count; i++ )
				var_choose_order[i] = i + 1;
		}
		else if ( schema_type == "bivium_Beginning2" ) {
			for ( unsigned i = 0; i < full_mask_var_count; i++ )
				var_choose_order[i] = i + 84 + 1;
		}
		else if ( schema_type == "bivium_Beginning_halved" ) {
			for ( unsigned i = 0; i < full_mask_var_count / 2; i++ )
				var_choose_order[k++] = i + 1;
			for ( unsigned i = 0; i < full_mask_var_count / 2; i++ )	
				var_choose_order[k++] = i + 84 + 1;
		}
		else if ( schema_type == "bivium_Ending1" ) {
			for ( unsigned i = 0; i < full_mask_var_count; i++ )
				var_choose_order[i] = 84 - i;
		}
		else if ( schema_type == "bivium_Ending2" ) {
			for ( unsigned i = 0; i < full_mask_var_count; i++ )
				var_choose_order[i] = 177 - i;
		}
		else if ( schema_type == "bivium_Ending_halved" ) {
				for ( unsigned i = 0; i < full_mask_var_count / 2; i++ )
					var_choose_order[k++] = 84 - i;
				for ( unsigned i = 0; i < full_mask_var_count / 2; i++ )	
					var_choose_order[k++] = 177 - i;
		}		
		else if ( schema_type == "a5" ) {
			for ( unsigned i = 0; i < 9; i++ )
				var_choose_order[k++] = i + 1;
			for ( unsigned i = 0; i < 11; i++ )
				var_choose_order[k++] = i + 19 + 1;
			for ( unsigned i = 0; i < 11; i++ )
				var_choose_order[k++] = i + 41 + 1;
		}
	}

	cout << "final var_choose_order.size() " << var_choose_order.size() << endl;
	sort( var_choose_order.begin(), var_choose_order.end() );
	cout << "var_choose_order" << endl;
	for ( unsigned i = 0; i < var_choose_order.size(); i++ )
		cout << var_choose_order[i] << " ";
	cout << endl;

	full_var_choose_order = var_choose_order;
	return true;
}

//---------------------------------------------------------
bool MPI_Base :: ReadVarCount( )
{
// Reading actual vatiable count CNF from file (skipping line "p cnf ...").
	unsigned int current_clause_count = 0,
				 current_lit_count = 0,
				 line_str_len = 0,
				 current_var_count = 0,
				 i = 0, k = 0,
				 lit_positiv_val = 0;
	string line_str, word_str;
	bool IsUncorrectLine = false;
	int lit_val, sign, lit_index, val;
	
	// check file with main CNF
	ifstream main_cnf( input_cnf_name, ios::in );
    if ( !main_cnf ) {
		std :: cerr << endl << "Error in opening of file with input CNF " 
			        << input_cnf_name << endl;
		return false;
	}

	// step 1 - get var_count and clause_count;
    while ( getline( main_cnf, line_str ) ) {
		if ( ( line_str[0] == 'p' ) || ( line_str[0] == 'c' ) )
			continue;

		else { // try to read line with clause 
			current_lit_count = 0;
			line_str = " " + line_str; // add space to line for correct work of parser
			for ( i = 0; i < line_str.length( ) - 1; i++ ) {
				IsUncorrectLine = false;
				if ( ( line_str[i] == ' ' ) && ( line_str[i + 1] != ' ' ) && 
					 ( line_str[i + 1] != '0' ) )
				{
					word_str = ""; // init string for cuttenr word from string
					k = i + 1; // k = index of first symbol in current word
					do {
						word_str += line_str[k];
						k++;
						if ( k == line_str.length( ) ) { // skip empty or uncorrect line
							/*std :: cout << "\n***In ReadVarCount skipped line " << 
								           line_str << endl;*/
							IsUncorrectLine = true;
							break;
						}
					} while ( line_str[k] != ' ' );

					if ( IsUncorrectLine )
						break;
					
					current_lit_count++;
					
					// convert value of literal to positiv int
					lit_positiv_val = abs( atoi( word_str.c_str( ) ) ); 
					if ( lit_positiv_val > current_var_count )
						current_var_count = lit_positiv_val;
				} // if ( ( line_str[i] == ' ' ) && ...
			} // for ( i = 0; i < line_str.length( ) - 1; i++ )

			if ( current_lit_count ) // if at least one lit exists then inc clause count
				current_clause_count++;
		}
	}

	var_count    = current_var_count;
	lit_count    = var_count * 2;
	clause_count = current_clause_count;

	clause_lengths.resize( clause_count );

	main_cnf.close( ); // reopen file
	main_cnf.clear( );
	main_cnf.open( input_cnf_name );
	current_clause_count = 0;

	// step 2 - get arrays of lengths
	while ( getline( main_cnf, line_str ) ) {
		if ( ( line_str[0] == 'p' ) || ( line_str[0] == 'c' ) )
			continue;

		else { // try to read line with clause 
			current_lit_count = 0;
			line_str = " " + line_str; // add space to line for correct work of parser
			for ( i = 0; i < line_str.length( ) - 1; i++ ) {
				IsUncorrectLine = false;
				if ( ( line_str[i] == ' ' ) && ( line_str[i + 1] != ' ' ) && 
					 ( line_str[i + 1] != '0' ) )
				{
					word_str = ""; // init string for cuttenr word from string
					k = i + 1; // k = index of first symbol in current word
					do {
						word_str += line_str[k];
						k++;
						if ( k == line_str.length( ) ) { // skip empty or uncorrect line
							/*std :: cout << "\n***In ReadVarCount skipped line " << 
								           line_str << endl;*/
							IsUncorrectLine = true;
							break;
						}
					} while ( line_str[k] != ' ' );

					if ( IsUncorrectLine )
						break;

					current_lit_count++;
				}
			} // for ( i = 0; i < line_str.length( ) - 1; i++ )

			if ( current_lit_count ) {
				current_clause_count++;
				clause_lengths[current_clause_count - 1] = current_lit_count;
			}
		}
	}

	main_cnf.close( );
	
	return true;
}


//---------------------------------------------------------
bool MPI_Base :: ReadIntCNF( )
{
// Reading CNF from file.
// Vars = {1, -1, 2, -2, ...} Literels = {2, 3, 4, 5, ...}
// clause_array        - array of clauses of literals
// lits_clause_array   - array of literals' clause indexes
// clause_num          - number of clauses
// var_count           - number of vars
// lits_num            - number of literals
// clause_lengths      - array of count of literals in clauses
// lits_clause_lengths - array of lengths of lits_clauses
    unsigned int current_clause_count = 0,
				 current_lit_count    = 0,
				 line_str_len         = 0,
				 current_var_count    = 0,
				 k                    = 0, 
				 val				  = 0,
				 sign				  = 0;
	unsigned input_var_num,
			 first_obj_var;
	int lit_val;
	string line_str, 
		   word_str;
	bool IncorrectLine;
	vector<int> :: iterator vec_it;

	if ( !ReadVarCount( ) ) {
		cerr << "Error in ReadVarCount" << endl;
		return false;
	}

	cout << "ReadVarCount() done" << endl;

	clause_array.resize( clause_count );

	int *lits_clause_current = new int[lit_count];

	for ( unsigned i = 0; i < clause_array.size(); ++i )
		clause_array[i].resize( clause_lengths[i] );

	// check file with main CNF
	ifstream main_cnf( input_cnf_name, ios::in );
    if ( !main_cnf ) {
		cerr << "Error in opening of file with input CNF with name" 
			 << input_cnf_name << endl;
		return false;
	}

	first_obj_var = 0;
	current_clause_count = 0;

	stringstream sstream;
	string str1, str2, str3, str4;
	int intval;
	bool Is_InpVar = false, Is_ConstrLen = false, Is_ObjLen = false, Is_ObjVars = false;
	while ( getline( main_cnf, line_str ) ) {
		if ( line_str[0] == 'p' )
			continue;
		if ( line_str[0] == 'c' ) { // in comment string can exist count of input variables
			//parse string for ex. "c 1452 input variables" or "c input variables 1452"
			sstream << line_str;
			sstream >> str1 >> str2;
			
			/*if ( str2 == "rslos" ) {
				while ( sstream >> intval )
					rslos_lengths.push_back( intval );
				cout << "rslos_count " << rslos_lengths.size() << endl;
				cout << "rslos lens: ";
				for ( unsigned i = 0; i < rslos_lengths.size(); i++ )
					cout << rslos_lengths[i] << " ";
				cout << endl;
				
				continue;
			}*/

			sstream >> str3 >> str4; // get and parse words in string
			sstream.str( "" ); sstream.clear( );

			if ( ( !Is_InpVar ) && ( str2 == "input" ) && ( str3 == "variables" ) ) {
				istringstream( str4 ) >> input_var_num;
				if ( input_var_num > 0 ) {
				    core_len = input_var_num;
				    if ( (core_len > MAX_CORE_LEN) || (core_len <= 0) ) {
						core_len = MAX_CORE_LEN;
						cout << "Warning. core_len > MAX_CORE_LEN or <= 0. Changed to MAX_CORE_LEN" << endl;
						cout << "core_len " << core_len << " MAX_CORE_LEN " << MAX_CORE_LEN << endl;
				    }
					for ( unsigned i=0; i < core_len; ++i )
						var_choose_order.push_back( i+1 );
				    Is_InpVar = true;
				    continue;
				}
			}
			if ( str2 == "var_set" ) {
				sstream << line_str;
				cout << "line_str " << line_str << endl;
				sstream >> str1; // remove "c"
				sstream >> str2; // remove "var_set"
				while ( sstream >> val ) {
					cout << val << " ";
					var_choose_order.push_back( val );
				}
				cout << endl;
				sstream.clear(); sstream.str();
				sort( var_choose_order.begin(), var_choose_order.end() );
				cout << "After reading var_set" << endl;
				cout << "var_choose_order.size() " << var_choose_order.size() << endl;
				for ( unsigned i=0; i < var_choose_order.size(); ++i )
					cout << var_choose_order[i] << " ";
				cout << endl;
				core_len = var_choose_order.size();
				cout << "core_len changed to " << core_len << endl;
			}
			
			if ( !Is_ConstrLen ) {
				if ( ( str2 == "constraint" ) && ( str3 == "clauses" ) ) {
					istringstream( str4 ) >> constr_clauses_count; 
					Is_ConstrLen = true;
					continue;
				}
			}
			if ( !Is_ObjLen ) {
				if ( ( ( str2 == "object" ) || ( str2 == "obj" ) ) && ( str3 == "clauses" ) ) {
					istringstream( str4 ) >> obj_clauses_count;
					Is_ObjLen = true;
					continue;
				}
			}
			if ( !Is_ObjVars ) {
				if ( ( str2 == "obj" ) && ( str3 == "vars" ) ) 
					istringstream( str4 ) >> first_obj_var;
				else if ( ( str2 == "object" ) && ( str3 == "variables" ) ) 
					istringstream( str4 ) >> first_obj_var;
				if ( first_obj_var > 0 ) {
					if ( !Is_ConstrLen ) constr_clauses_count = clause_count;

					// read obj vars indexes from string "c obj vars ...")
					obj_vars_count = 0;
					line_str += ' ';
					while ( line_str.length( ) ) {
						while ( ( line_str.length( ) > 0 ) && ( line_str[0] == ' ' ) ) 
							line_str.erase( 0, 1 ); // skip spaces
						if ( !line_str.length( ) )
							break;
						int space_pos = line_str.find( " " ); // find space after word
						word_str = line_str.substr( 0, space_pos ); // get word
						line_str.erase( 0, space_pos + 1 ); // delete word and space after it
						int word_value = atoi( word_str.c_str( ) );
						if ( word_value > 0 ) { // if it is number > 0
							if ( obj_vars_count < 32 ) {
								obj_vars[obj_vars_count] = word_value; // add to array
								obj_vars_count++;
							}
							else
								cerr << "String c obj vars ... contains too many values ";
						}
					}
					Is_ObjVars = true;
					continue;
				} 
			}  // if ( !Is_ObjVars )
		}
		else // if ( ( line_str[0] != 'p' ) && ( line_str[0] != 'c' ) )
		{
			// try to read line with clause
			current_lit_count = 0; // cuttenr count of lits in current clause
			line_str = " " + line_str;
			for ( unsigned i = 0; i < line_str.length( ) - 1; i++ ) {
				IncorrectLine = false;
				if ( ( line_str[i] == ' ' ) && ( line_str[i + 1] != ' ' ) && 
					 ( line_str[i + 1] != '0' ) )
				{
					word_str = "";
					k = i + 1;
					do {
						word_str += line_str[k];
						k++;
						if ( k == line_str.length( ) ) { // skip empty or uncorrect line
							/*std :: cout << "\n***In ReadVarCount skipped line " << 
								           line_str << endl;*/
							IncorrectLine = true;
							break;
						}
					} while ( line_str[k] != ' ' );

					if ( IncorrectLine )
						break;

					lit_val = atoi( word_str.c_str( ) ); // word -> lit value
					if ( !lit_val ) { // if non-number or '0' (lit > 0) then rerurn error;
						cout << "\n Error in ReadIntCNF. literal " << word_str << " is non-number";
						return false;
					}
					else if ( lit_val < 0 ) {
						lit_val = -lit_val; 
						sign = 1;
					}
					else if ( lit_val > 0 )
						sign = 0;

					// fill attay of clauses
					val = ( lit_val << 1 ) + sign; // literal value, 1 -> 2, -1 -> 3, 2 -> 4, -2 -> 5
					clause_array[current_clause_count][current_lit_count] = val;
					current_lit_count++;
				} // if ( ( line_str[i] == ' ' ) ...
			} // for ( i = 0; i < line_str.length( ) - 1; i++ )
			
			if ( current_lit_count )
				current_clause_count++;
		} 
	} // while ( getline( main_cnf, line_str ) )
	
	// if PB data is correct then turn on PB mode
	/*if ( ( constr_clauses_count > 0 ) && ( obj_vars_count > 0 ) ) 
	{
		IsPB = true;
		solver_type = 3;
		if ( ( best_lower_bound > -1 ) && ( upper_bound > 0 ) ) // if bounds then equality mode
		{
			PB_mode = 2;
			if ( best_lower_bound > upper_bound )  // check correctness
				best_lower_bound = upper_bound;
		}
		else PB_mode = 1;
	}*/

	main_cnf.close( );
	
	full_var_choose_order = var_choose_order;
	// fill indexes of core variables
	k=0;
	for ( vec_it = full_var_choose_order.begin(); vec_it != full_var_choose_order.end(); ++vec_it )
		core_var_indexes.insert(pair<int,unsigned>( *vec_it, k++ ));
	if ( verbosity > 0 ) {
		cout << "core_var_indexes" << endl;
		for ( map<int,unsigned> :: iterator map_it = core_var_indexes.begin(); map_it != core_var_indexes.end(); ++map_it )
			cout << map_it->first << " " << map_it->second << endl;
	}
	
	cout << "ReadIntCNF() done" << endl;
	
	return true;
}

//---------------------------------------------------------
bool MPI_Base :: CheckSATset( vector<int> &lit_SAT_set_array )
{
// Check given SAT set
	int cnf_sat_val = 1,
		i = -1,
		checked_clauses = 0,
		j = -1,
		current_lit_val = -1,
		var_index = -1,
		current_sat_val = -1,
		clause_sat_val = -1;
	//
	i = 0;
	// check clauses while all of checked is SAT
	while ( ( cnf_sat_val ) && ( checked_clauses != clause_count ) ) {
		// add one to count of checked clauses in main CNF
        checked_clauses++;
		clause_sat_val = 0;
		j = 0;
		// check literals while all of checked is UNSAT
		while ( ( !clause_sat_val ) && ( j < clause_lengths[i] ) ) {
			current_lit_val = clause_array[i][j];
			var_index = current_lit_val/2;
			current_sat_val = lit_SAT_set_array[var_index - 1]; 
			if ( current_lit_val == current_sat_val )
				clause_sat_val++;
			j++;
		} // while ( ( !clause_sat_val ) && 
		//
		// if we found one UNSAT clause CNF is UNSAT too
		if ( !clause_sat_val )
			--cnf_sat_val;
		i++;
	} // while ( ( cnf_sat_val ) && ( checked_clauses != clause_num ) )
	return ( cnf_sat_val > 0 ) ? true : false;
}

//---------------------------------------------------------
bool MPI_Base :: AnalyzeSATset( )
{
// Reading of SAT set and check it for correction
	int lits_num = 0,
		int_answer = 0,
		sign = 0,
		val = 0,
		str_output_cur_ind = 0;
	unsigned answer_var_count;
	bool bIsSATSetExist = false;
	string answer_file_name,
		   output_file_name = "output",
		   str_answer,
		   line_buffer;
	ofstream answer_file,
			 output_file;
	vector<int> lit_SAT_set_array;
	stringstream sstream;
	
	cout << "Start of AnalyzeSATset" << endl;
	lit_SAT_set_array.resize( b_SAT_set_array.size() );
	for ( unsigned i = 0; i < lit_SAT_set_array.size(); ++i )
		lit_SAT_set_array[i] = ( ( i + 1 ) << 1 ) + (b_SAT_set_array[i] ? 0 : 1);
	// if SAT set exist then check it
	if ( verbosity > 0 )
		cout << "Before CheckSATset()" << endl;
	if ( !CheckSATset( lit_SAT_set_array ) ) {
		cerr << "Error in checking of SAT set" << endl;
		return false;
	}
	if ( verbosity > 0 )
		cout << "CheckSATset() done" << endl;
	
	// open file for writing of answer in zchaff style
	answer_file_name = "sat_sets_";
	answer_file_name += input_cnf_name;
	
	unsigned k = 1;
	for( unsigned i = 1; i < answer_file_name.length( ); ++i ) {
		if ( ( answer_file_name[i] == '.' ) || ( answer_file_name[i] == '\\' ) || 
			 ( answer_file_name[i] == '/' ) || ( answer_file_name[i] == '//' ) ||
			 ( answer_file_name[i] == ':' ) || 
			 ( ( answer_file_name[i] == '_' ) && ( answer_file_name[i - 1] == '_' ) ) )
			continue;
		else {
			answer_file_name[k] = answer_file_name[i];
			k++;
		}
	}
	answer_file_name.resize( k );
	sstream << "_" << rank;
	answer_file_name += sstream.str();
	sstream.clear(); sstream.str("");

	cout << "answer_file_name " << answer_file_name << endl; 
	answer_file.open( answer_file_name.c_str( ), ios :: app );
	if ( !( answer_file.is_open( ) ) ) {
		std :: cerr << "Error in opening of file with answer in zchaff style " 
			        << answer_file_name << endl;
		return false;
	}
	
	answer_var_count = core_len;
	
	sstream << "SAT" << endl;
	for ( unsigned i = 0; i < b_SAT_set_array.size(); ++i )
		sstream << b_SAT_set_array[i];
	sstream << endl; 
	answer_file << sstream.rdbuf( );
	answer_file.close( );
	lit_SAT_set_array.clear();
	
	return true;
}

// TODO move to Addit_func namespace
unsigned MPI_Base ::  uint_rand() {
	boost::random::uniform_int_distribution<uint32_t> dist;
	return dist(gen);
}

bool MPI_Base ::  bool_rand() {
	boost::random::uniform_int_distribution<uint32_t> dist(0,1);
	return bool( dist(gen) );
}

void MPI_Base :: MakeRandArr( vector< vector<unsigned> > &rand_arr, unsigned vec_len, unsigned rnd_uint32_count )
{
// make array of pseudorandom values using Mersenne Twister generator
	rand_arr.resize( vec_len );
	vector< vector<unsigned> > :: iterator it;
	for ( it = rand_arr.begin(); it != rand_arr.end(); it++ ) {
		(*it).resize( rnd_uint32_count );
		for ( unsigned j = 0; j < (*it).size(); j++ )
			(*it)[j] = uint_rand();
	}
}

void MPI_Base :: MakeUniqueRandArr( vector<unsigned> &rand_arr, unsigned rand_arr_len, 
							        unsigned max_rand_val )
{
// make array of different pseudorandom values
	if ( max_rand_val < rand_arr_len )
		max_rand_val = rand_arr_len;
	unsigned rand_numb;
	rand_arr.resize( rand_arr_len );
	
	bool IsOldValue;
	for ( unsigned i = 0; i < rand_arr_len; i++ ) {
		do { // if value is not unique get value again
			rand_numb = uint_rand();
			rand_numb %= max_rand_val;
			IsOldValue = false;
			for ( unsigned k = 0; k < i; ++k ) {
				if ( rand_numb == rand_arr[k] ) {
					IsOldValue = true;
					break;
				}
			}
		} while ( IsOldValue );
		rand_arr[i] = rand_numb; // new values
	}
}