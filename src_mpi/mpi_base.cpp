#include "mpi_base.h"

#pragma warning( disable : 4996 )

const int    MAX_WORD_LENGTH			   = 64;
const int    MAX_LINE_LENGTH               = 524288;
const int    MAX_LINE_LENGTH_2             = 8192;
const int    MEDIUM_STRING_LEN             = 4096;
const double TRANSP_COAST                  = 0.000001;
const int    NUM_KEY_BITS                  = 64;

using namespace std;

//=================================================================================================
// Constructor/Destructor:

MPI_Base :: MPI_Base( ) :
	sat_count            ( 0 ),
	corecount			 ( 0 ),
	solver_name          ( "" ),
	core_len             ( 0 ),
	koef_val             ( 4 ),
	schema_type          ( "" ),
	var_count            ( 0 ),
	clause_count         ( 0 ),
	full_mask_var_count  ( 0 ),
	start_activity		 ( 0 ),
	part_mask_var_count  ( 0 ),
	all_tasks_count		 ( 0 ),
	isConseq             ( false ),
	cnf_in_set_count     ( 100 ),
	verbosity			 ( 0 ),
	known_point_file_name ( "known_point" ),
	isSolveAll           ( false ),
	isPredict            ( false ),
	max_solving_time     ( 0 ),
	max_nof_restarts     ( 0 ),
	first_stream_var_index ( 0 ),
	te ( 0 ),
	known_bits ( 0 ),
	output_len ( 0 ),
	isMakeSatSampleAnyWay ( false ),
	isSolverSystemCalling ( false ),
	process_sat_count ( 0 ),
	known_vars_count ( 0 ),
	isPlainText (false),
	evaluation_type("time"),
	rank (0),
	total_start_time (0.0),
	base_known_vars_count(0)
{
	for ( unsigned i = 0; i < FULL_MASK_LEN; i++ )
		full_mask[i] = part_mask[i] = mask_value[i] = 0;
}

MPI_Base :: ~MPI_Base( )
{ }

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
		full_mask[cur_uint_ind] += 1 << ( var_index % UINT_LEN );
	}

	// get first full_mask_var_count vars from  array var_choose_order
	for ( unsigned i = 0; i < part_mask_var_count; ++i ) {
		var_index = var_choose_order[i] - 1;
		cur_uint_ind = var_index / UINT_LEN;
		part_mask[cur_uint_ind] += 1 << ( var_index % UINT_LEN );
	}
	
	if ( verbosity > 1 ) {
		cout << "made part_mask" << endl;
		for( unsigned i = 0; i < FULL_MASK_LEN; ++i )
			cout << part_mask[i] << " ";
		cout << endl;
		unsigned bit_count = 0;
		for ( unsigned j=0; j < FULL_MASK_LEN; ++j )
			bit_count += BitCount( part_mask[j] );
		cout << "part_mask bit_count " << bit_count << endl;
	}
	
	if ( verbosity > 1 )
		cout << "GetMainMasksFromVarChoose() end" << endl;
	
	return true;
}

bool MPI_Base :: MakeAssignsFromMasks( unsigned *full_mask, 
									   unsigned *part_mask, 
									   unsigned *mask_value,
									   vec<Lit> &known_dummy,
									   vec< vec<Lit> > &dummy_vec,
									   int rank )
{
    // convert masks to a vector of Literals
	unsigned variate_vars_count = 0, known_vars_count = 0;
	full_mask_var_count = 0;
	for ( unsigned i = 0; i < FULL_MASK_LEN; i++ ) {
		variate_vars_count  += BitCount( full_mask[i] ^ part_mask[i] );
		full_mask_var_count += BitCount( full_mask[i] );
	}
	known_vars_count = full_mask_var_count - variate_vars_count;
	
	unsigned mask;
	int cur_var_ind;
	
	// get literals which can be assigned once (from part_mask)
	for (unsigned i = 0; i < FULL_MASK_LEN; i++) {
		for (unsigned j = 0; j < UINT_LEN; j++) {
			mask = (1 << j);
			cur_var_ind = i * UINT_LEN + j;
			if (part_mask[i] & mask) // one common vector send by control process
				known_dummy.push((mask_value[i] & mask) ? mkLit(cur_var_ind) : ~mkLit(cur_var_ind) );
		}
	}

	if (known_dummy.size() != known_vars_count) {
		cerr << "known_dummy.size() != known_vars_count" << endl;
		cerr << known_dummy.size() << " != " << known_vars_count << endl;
#ifdef _MPI
		MPI_Abort( MPI_COMM_WORLD, 0 );
#endif
	}
	
	if (verbosity > 2)
		cout << "known_dummy.size() " << known_dummy.size() << endl;

	// determine the number of assumptions and their lengths
	unsigned problems_count = 1 << variate_vars_count;
	dummy_vec.resize(problems_count);
	for (int i = 0; i < dummy_vec.size(); i++)
		dummy_vec[i].resize(variate_vars_count);

	// get array of literals which must be checked
	for ( unsigned lint = 0; lint < problems_count; lint++ ) {
		unsigned range_mask;
		unsigned range_mask_ind = 0;
		unsigned index = 0;
		for ( unsigned i = 0; i < FULL_MASK_LEN; i++ ) {
			for ( unsigned j = 0; j < UINT_LEN; j++ ) {
				mask = ( 1 << j );
				cur_var_ind = i * UINT_LEN + j;
				if ( ( full_mask[i] & mask ) && !( part_mask[i] & mask ) ) { // set of vectors
					range_mask = 1 << range_mask_ind;
					range_mask_ind++;
					dummy_vec[lint][index++] = (lint & range_mask) ? mkLit(cur_var_ind) : ~mkLit(cur_var_ind);
				}
			}
		}
	}

	if (verbosity > 2)
		cout << "dummy_vec.size() " << dummy_vec.size() << endl;

	if (verbosity > 0)
		cout << "MakeAssignsFromMasks() done " << endl;

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
		for( unsigned i = 0; i < FULL_MASK_LEN; i++ ) {
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

	return true;
}

//---------------------------------------------------------
bool MPI_Base::getValuesFromIntegers(vector<vector<int>> cartesian_elements)
{
	unsigned mask;
	int cur_integer_index;
	unsigned cur_var_index = 0;
	unsigned k;
	
	for (unsigned cartesian_element_index = 0; cartesian_element_index < cartesian_elements.size(); cartesian_element_index++) {
		k = 0;
		cur_integer_index = cartesian_elements[cartesian_element_index][k];
		for (unsigned i = 0; i < FULL_MASK_LEN; i++) {
			if (k == cartesian_elements[cartesian_element_index].size())
				break;
			for (unsigned j = 0; j < UINT_LEN; j++) {
				cur_var_index = i*UINT_LEN + j;
				mask = (1 << j);
				if (cur_var_index == cur_integer_index) {
					values_arr[cartesian_element_index][i] += mask;
					k++;
					if (k == cartesian_elements[cartesian_element_index].size())
						break;
					cur_integer_index = cartesian_elements[cartesian_element_index][k];
				}
			} // for( j = 0; j < UINT_LEN; j++ )
		}
	} // for( lint = 0; lint < part_var_power; lint++ )

	return true;
}

//---------------------------------------------------------
bool MPI_Base :: makeStandardMasks( unsigned &part_var_power )
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
bool MPI_Base::makeIntegerMasks(vector<vector<int>> cartesian_elements)
{
	if (!GetMainMasksFromVarChoose(var_choose_order)) {
		cout << "Error in GetMainMasksFromVarChoose" << endl;
		return false;
	}
	cout << "Correct end of GetMasksFromVarChoose" << endl;

	if (!getValuesFromIntegers(cartesian_elements)) {
		cout << "Error in getValuesFromIntegers" << endl;
		return false;
	}
	
	cout << "Correct end of makeIntegerMasks" << endl;
	
	return true;
}

//---------------------------------------------------------
bool MPI_Base :: MakeVarChoose()
{
// Make array var_choose_order with vars sorted by given rule
	string str;
	stringstream sstream;
	int val;

	cout << "start MakeVarChoose()" << endl;
	
	// if file with decomposition set exists
	ifstream known_point_file( known_point_file_name.c_str() );
	if ( known_point_file.is_open() ) {
		cout << "reading known_point_file " << known_point_file_name << endl;
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
	else {
		var_choose_order.resize(full_var_choose_order.size());
		int k = 0;
		if ( (schema_type == "0") || (schema_type == "") ) { // == bivium_Beginning1
			for (unsigned i = 0; i < full_var_choose_order.size(); i++)
				var_choose_order[i] = full_var_choose_order[i];
		}
		else if (schema_type == "bivium_Beginning2") {
			for (unsigned i = 0; i < full_mask_var_count; i++)
				var_choose_order[i] = i + 84 + 1;
		}
		else if (schema_type == "bivium_Beginning_halved") {
			for (unsigned i = 0; i < full_mask_var_count / 2; i++)
				var_choose_order[k++] = i + 1;
			for (unsigned i = 0; i < full_mask_var_count / 2; i++)
				var_choose_order[k++] = i + 84 + 1;
		}
		else if (schema_type == "bivium_Ending1") {
			for (unsigned i = 0; i < full_mask_var_count; i++)
				var_choose_order[i] = 84 - i;
		}
		else if (schema_type == "bivium_Ending2") {
			for (unsigned i = 0; i < full_mask_var_count; i++)
				var_choose_order[i] = 177 - i;
		}
		else if (schema_type == "bivium_Ending_halved") {
			for (unsigned i = 0; i < full_mask_var_count / 2; i++)
				var_choose_order[k++] = 84 - i;
			for (unsigned i = 0; i < full_mask_var_count / 2; i++)
				var_choose_order[k++] = 177 - i;
		}
		else if (schema_type == "a5") {
			for (unsigned i = 0; i < 9; i++)
				var_choose_order[k++] = i + 1;
			for (unsigned i = 0; i < 11; i++)
				var_choose_order[k++] = i + 19 + 1;
			for (unsigned i = 0; i < 11; i++)
				var_choose_order[k++] = i + 41 + 1;
		}
	}

	cout << "var_choose_order.size() " << var_choose_order.size() << endl;
	sort( var_choose_order.begin(), var_choose_order.end() );
	cout << "var_choose_order" << endl;
	for ( unsigned i = 0; i < var_choose_order.size(); i++ )
		cout << var_choose_order[i] << " ";
	cout << endl;

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
	
	// check file with main CNF
	ifstream main_cnf( input_cnf_name.c_str(), ios::in );
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
	main_cnf.open( input_cnf_name.c_str() );
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
bool MPI_Base::ReadIntCNF()
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
		current_lit_count = 0,
		line_str_len = 0,
		current_var_count = 0,
		k = 0,
		val = 0,
		sign = 0;
	unsigned first_obj_var;
	int lit_val;
	string line_str,
		word_str;
	bool IncorrectLine;
	vector<int> ::iterator vec_it;

	cout << "Start of ReadIntCNF()" << endl;
	if (!ReadVarCount()) {
		cerr << "Error in ReadVarCount" << endl; return false;
	}

	cout << "ReadVarCount() done" << endl;

	clause_array.resize(clause_count);

	for (unsigned i = 0; i < clause_array.size(); ++i)
		clause_array[i].resize(clause_lengths[i]);

	// check file with main CNF
	ifstream main_cnf(input_cnf_name.c_str(), ios::in);
	if (!main_cnf.is_open()) {
		cerr << "Error in opening of file with input CNF with name"
			<< input_cnf_name << endl;
		exit(1);
	}

	first_obj_var = 0;
	current_clause_count = 0;

	stringstream sstream;
	string str1, str2, str3, str4, str5;
	unsigned ui;
	bool Is_InpVar = false, Is_ConstrLen = false, Is_ObjLen = false, Is_ObjVars = false;
	while (getline(main_cnf, line_str)) {
		if (line_str[0] == 'p')
			continue;
		if (line_str[0] == 'c') { // in comment string can exist count of input variables
			//parse string for ex. "c 1452 input variables" or "c input variables 1452"
			sstream.str(""); sstream.clear();
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

			if (str2 == "known_bits") {
				istringstream(str3) >> known_bits;
				cout << "known_bits " << known_bits << endl;
				continue;
			}
			if ((str2 == "output") && (str3 == "variables")) {
				istringstream(str4) >> output_len;
				continue;
			}
			if ((str2 == "output") && (str3 == "vars")) {
				sstream >> str5;
				istringstream(str5) >> output_len;
				continue;
			}
			if (!Is_InpVar) {
				if ((str2 == "input") && (str3 == "variables"))
					istringstream(str4) >> ui;
				else if ((str2 == "input") && (str3 == "vars")) {
					sstream >> str5;
					istringstream(str5) >> ui;
				}
				cout << "input variables " << ui << endl;
				if (!core_len)
					core_len = ui; // if core_len didn't set manually, read from file
				if ((core_len > MAX_CORE_LEN) || (core_len <= 0)) {
					core_len = MAX_CORE_LEN;
					cout << "Warning. core_len > MAX_CORE_LEN or <= 0. Changed to MAX_CORE_LEN" << endl;
					cout << "core_len " << core_len << " MAX_CORE_LEN " << MAX_CORE_LEN << endl;
				}
				Is_InpVar = true;
				continue;
			}
			sstream.str(""); sstream.clear();

			if (str2 == "var_set") {
				sstream << line_str;
				cout << "line_str " << line_str << endl;
				sstream >> str1; // remove "c"
				sstream >> str2; // remove "var_set"
				while (sstream >> val) {
					cout << val << " ";
					full_var_choose_order.push_back(val);
				}
				cout << endl;
				sstream.clear(); sstream.str();
				sort(full_var_choose_order.begin(), full_var_choose_order.end());
				cout << "After reading var_set" << endl;
				cout << "full_var_choose_order.size() " << full_var_choose_order.size() << endl;
				for (unsigned i = 0; i < full_var_choose_order.size(); ++i)
					cout << full_var_choose_order[i] << " ";
				cout << endl;
				core_len = full_var_choose_order.size();
				cout << "core_len changed to " << core_len << endl;
			}
		}
		else // if ( ( line_str[0] != 'p' ) && ( line_str[0] != 'c' ) )
		{
			// try to read line with clause
			current_lit_count = 0; // current count of lits in current clause
			line_str = " " + line_str;
			for (unsigned i = 0; i < line_str.length() - 1; i++) {
				IncorrectLine = false;
				if ((line_str[i] == ' ') && (line_str[i + 1] != ' ') &&
					(line_str[i + 1] != '0'))
				{
					word_str = "";
					k = i + 1;
					do {
						word_str += line_str[k];
						k++;
						if (k == line_str.length()) { // skip empty or uncorrect line
							/*std :: cout << "\n***In ReadVarCount skipped line " <<
										   line_str << endl;*/
							IncorrectLine = true;
							break;
						}
					} while (line_str[k] != ' ');

					if (IncorrectLine)
						break;

					lit_val = atoi(word_str.c_str()); // word -> lit value
					if (!lit_val) { // if non-number or '0' (lit > 0) then rerurn error;
						cout << "\n Error in ReadIntCNF. literal " << word_str << " is non-number";
						return false;
					}
					else if (lit_val < 0) {
						lit_val = -lit_val;
						sign = 1;
					}
					else if (lit_val > 0)
						sign = 0;

					// fill attay of clauses
					val = (lit_val << 1) + sign; // literal value, 1 -> 2, -1 -> 3, 2 -> 4, -2 -> 5
					clause_array[current_clause_count][current_lit_count] = val;
					current_lit_count++;
				} // if ( ( line_str[i] == ' ' ) ...
			} // for ( i = 0; i < line_str.length( ) - 1; i++ )

			if ((te > 0) && (current_lit_count == 1))
				cout << "Warning. ( te > 0 ) && ( current_lit_count == 1 ). change CNF file to template one" << endl;

			if (current_lit_count == 1)
				known_vars_count++;

			if (current_lit_count)
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

	main_cnf.close();

	// fill indexes of core variables
	/*k=0;
	for ( vec_it = full_var_choose_order.begin(); vec_it != full_var_choose_order.end(); ++vec_it )
		core_var_indexes.insert(pair<int,unsigned>( *vec_it, k++ ));
	if ( verbosity > 0 ) {
		cout << "core_var_indexes" << endl;
		for ( map<int,unsigned> :: iterator map_it = core_var_indexes.begin(); map_it != core_var_indexes.end(); ++map_it )
			cout << map_it->first << " " << map_it->second << endl;
	}*/

	if (core_len)
		nonoutput_len = core_len;
	else {
		nonoutput_len = var_count - output_len;
	}
	
	// if wasn't defined by var_set
	if (full_var_choose_order.empty()) {
		if (known_bits) {
			core_len -= known_bits;
			for (unsigned i = 0; i < core_len; i++)
				full_var_choose_order.push_back(i + 1);
			cout << "new core_len (less due to nonzero value of known_bits) " << core_len << endl;
		}
		else {
			for (unsigned i = 0; i < nonoutput_len; i++)
				full_var_choose_order.push_back(i + 1);
		}
	}
	else {
		nonoutput_len = core_len = full_var_choose_order.size();
	}
	
	if ( ( isPredict ) && ( !core_len ) ) {
		cerr << "core_len == 0 in predict mode" << endl;
		return false;
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
bool MPI_Base::AnalyzeSATset(double cnf_time_from_node, 
							 int dummy_index, 
							 vector<bool> &dummy_bool,
							 unsigned *full_mask,
							 unsigned *part_mask,
							 unsigned *value)
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
	lit_SAT_set_array.resize(b_SAT_set_array.size());
	for (unsigned i = 0; i < lit_SAT_set_array.size(); ++i)
		lit_SAT_set_array[i] = ((i + 1) << 1) + (b_SAT_set_array[i] ? 0 : 1);
	// if SAT set exist then check it
	if (verbosity > 0)
		cout << "Before CheckSATset()" << endl;
	if (!CheckSATset(lit_SAT_set_array)) {
		cerr << "Error in checking of SAT set" << endl;
		return false;
	}
	if (verbosity > 0)
		cout << "CheckSATset() done" << endl;

	// open file for writing of answer in zchaff style
	answer_file_name = "sat_sets_";
	answer_file_name += input_cnf_name;

	unsigned k = 1;
	for (unsigned i = 1; i < answer_file_name.length(); ++i) {
		if ((answer_file_name[i] == '.') || (answer_file_name[i] == '\\') ||
			(answer_file_name[i] == '/') || (answer_file_name[i] == ':') ||
			((answer_file_name[i] == '_') && (answer_file_name[i - 1] == '_')))
			continue;
		else {
			answer_file_name[k] = answer_file_name[i];
			k++;
		}
	}
	answer_file_name.resize(k);
	sstream << "_" << rank;
	answer_file_name += sstream.str();
	sstream.clear(); sstream.str("");

	cout << "answer_file_name " << answer_file_name << endl;
	answer_file.open(answer_file_name.c_str(), ios::app);
	if (!(answer_file.is_open())) {
		cerr << "Error in opening of file with answer in zchaff style "
			<< answer_file_name << endl;
		return false;
	}

	answer_var_count = core_len;

#ifdef _MPI
	sstream << "total_time_from_start " << MPI_Wtime() - total_start_time << " s" << endl;
	sstream << "dummy_index " << dummy_index << endl;
	if (dummy_bool.size() > 0) {
		sstream << "dummy " << dummy_index << endl;
		for (unsigned j=0; j < dummy_bool.size(); j++)
			sstream << (dummy_bool[j] ? 1 : 0);
		sstream << endl;
		for (unsigned t = 0; t < FULL_MASK_LEN; t++)
			if ((full_mask[t]) || (t == 0))
				sstream << full_mask[t] << " ";
		sstream << endl;
		for (unsigned t = 0; t < FULL_MASK_LEN; t++)
			if ((part_mask[t]) || (t == 0))
				sstream << part_mask[t] << " ";
		sstream << endl;
		for (unsigned t = 0; t < FULL_MASK_LEN; t++)
			if ((mask_value[t]) || (t == 0))
				sstream << mask_value[t] << " ";
		sstream << endl;
	}
	sstream << cnf_time_from_node << " s SAT " << endl;
	sstream << endl;
#endif
	for ( unsigned i = 0; i < b_SAT_set_array.size(); ++i )
		sstream << b_SAT_set_array[i];
	sstream << endl; 
	answer_file << sstream.rdbuf( );
	answer_file.close( );
	lit_SAT_set_array.clear();
	
	return true;
}

void MPI_Base :: MakeRandArr( vector< vector<unsigned> > &rand_arr, unsigned vec_len, unsigned rnd_uint32_count )
{
// make array of pseudorandom values using Mersenne Twister generator
	rand_arr.resize( vec_len );
	vector< vector<unsigned> > :: iterator it;
	for ( it = rand_arr.begin(); it != rand_arr.end(); it++ ) {
		(*it).resize( rnd_uint32_count );
		for ( unsigned j = 0; j < (*it).size(); j++ )
			(*it)[j] = uint_rand( gen );
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
			rand_numb = uint_rand( gen );
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

void MPI_Base::MakeSingleSatSample(
		vector<bool> &state_vec, 
		vector<bool> &stream_vec, 
		const int seed_num,
		Solver *S,
		vector<lbool> predefined_vars)
{
	boost::random::mt19937 gen_local;
	//gen_local.seed( static_cast<unsigned>(seed_num^1234567));
	gen_local.seed( static_cast<unsigned>(seed_num));
	for (int i=0;i<100000;i++) gen_local;

	if (predefined_vars.empty())
		predefined_vars.resize(core_len, l_Undef);
	assert( predefined_vars.size()==core_len);

	vec<Lit> dummy;
	for ( unsigned i=0; i < core_len; i++ ){
		state_vec.push_back(l_Undef==predefined_vars[i]? bool_rand(gen_local) : l_True==predefined_vars[i]);
		//state_vec.push_back(bool_rand(gen_local));
		dummy.push(~mkLit(i, state_vec[i]));
	}
	Minisat::lbool ret = S->solveLimited( dummy );

	if ( ret != l_True ) {
		cerr << "in makeSatSample() ret != l_True" << endl;
		exit(1);
	}

	for( int i=state_vec.size(); i < S->model.size() - (int)output_len; i++ )
		state_vec.push_back( (S->model[i] == l_True) ? true : false );
	for( int i=S->model.size() - output_len; i < S->model.size(); i++ )
		stream_vec.push_back( (S->model[i] == l_True) ? true : false );
}

void MPI_Base::MakeSatSample(vector< vector<bool> > &state_vec_vec,
							 vector< vector<bool> > &stream_vec_vec,
							 vector< vector<bool> > &plain_text_vec_vec,
							 int rank)
{
	stringstream sstream;
	sstream << "known_sat_sample_" << rank;
	string known_sat_sample_file_name = sstream.str();
	sstream.str(""); sstream.clear();
	fstream known_sat_sample_file(known_sat_sample_file_name, ios_base::in );
	vector<bool> state_vec, stream_vec, plain_text_vec;
	
	if ( ( isMakeSatSampleAnyWay ) || (!known_sat_sample_file.is_open()) ) { // empty file
	//if ( file.peek() == fstream::traits_type::eof() ) { // if file is empty
		// make [sample_size] different pairs <register_state, keystream> via generating secret keys
		cout << "file known_sat_sample is empty. making SAT sample" << endl;
		
		boost::random::mt19937 gen_known_vars;
		const unsigned known_vars_seed = static_cast <unsigned> (time(NULL)) + (unsigned)rank * 1000000;
		cout << "known_vars_seed " << known_vars_seed << endl;
		gen_known_vars.seed(known_vars_seed);
		
		if (state_vec_vec.size() == 0) {
			state_vec.resize(core_len);
			for (unsigned i = 0; i < cnf_in_set_count; i++) {
				for (unsigned j = 0; j < core_len; j++)
					state_vec[j] = bool_rand(gen_known_vars);
				state_vec_vec.push_back(state_vec);
			}
		}
				
		// additionally plaintext is nedded 
		if (isPlainText) {
			unsigned ciphertext_len = 0;
			if (input_cnf_name.find("32bits") != string::npos)
				ciphertext_len = 32;
			else if (input_cnf_name.find("64bits") != string::npos)
				ciphertext_len = 64;
			else if (input_cnf_name.find("128bits") != string::npos)
				ciphertext_len = 128;
			else if (input_cnf_name.find("256bits") != string::npos)
				ciphertext_len = 256;
			else if (input_cnf_name.find("512bits") != string::npos)
				ciphertext_len = 512;
			else {
				cerr << "ciphertext_len == 0" << endl;
				exit(1);
			}
			cout << "ciphertext_len " << ciphertext_len << endl;
			plain_text_vec.resize(ciphertext_len);
			for (unsigned i = 0; i < cnf_in_set_count; i++) {
				for (unsigned j = 0; j < ciphertext_len; j++)
					plain_text_vec[j] = bool_rand(gen);
				plain_text_vec_vec.push_back(plain_text_vec);
			}
			
			cout << "plain_text_vec_vec.size() " << plain_text_vec_vec.size() << endl;
			cout << "plain_text_vec_vec[0].size() " << plain_text_vec_vec[0].size() << endl;
			
			for (unsigned i = 0; i < plain_text_vec_vec.size(); i++)
				for (unsigned j = 0; j < plain_text_vec_vec[i].size(); j++)
					state_vec_vec[i].push_back(plain_text_vec_vec[i][j]);

			cout << "state_vec_vec.size() " << state_vec_vec.size() << endl;
			cout << "state_vec_vec[0].size() " << state_vec_vec[0].size() << endl;
		}
		
		// get state of additional variables
		Problem cnf;
		Solver *S;
		Minisat::lbool ret;
		minisat22_wrapper m22_wrapper;
		ifstream in( input_cnf_name.c_str() );
		m22_wrapper.parse_DIMACS_to_problem(in, cnf);
		in.close();
		S = new Solver();
		S->addProblem(cnf);
		vec<Lit> dummy;
		int cur_var_ind;

		//int state_vec_len = state_vec_vec[0].size();
		//for ( vector< vector<bool> > :: iterator x = state_vec_vec.begin(); x != state_vec_vec.end(); x++ ) 
		unsigned long long unsat_genereted_count = 0, undef_genereted_count = 0;
		for (auto &cur_state : state_vec_vec) {
			cur_var_ind = 0;
			for (unsigned i = 0; i < cur_state.size(); i++) {
				dummy.push(cur_state[i] ? mkLit( cur_var_ind ) : ~mkLit( cur_var_ind ) );
				cur_var_ind++;
			}
			ret = S->solveLimited( dummy );
			dummy.clear();
			if ( ret == l_True ) {
				for (int i = S->model.size() - output_len; i < S->model.size(); i++)
					stream_vec.push_back((S->model[i] == l_True) ? true : false);
				stream_vec_vec.push_back(stream_vec);
				stream_vec.clear();
			}
			else if (ret == l_False)
				unsat_genereted_count++;
			else if (ret == l_Undef)
				undef_genereted_count++;
			if ((unsat_genereted_count) && (unsat_genereted_count % 100000 == 0)) {
				cout << "unsat_generated_count " << unsat_genereted_count << endl;
				//cout << "undef_generated_count " << undef_genereted_count << endl;
				cout << "state_vec_vec.size() " << state_vec_vec.size() << endl;
			}
		}
	
		cout << "unsat_generated_count " << unsat_genereted_count << endl;
		sstream << "state" << endl;
		for ( vector< vector<bool> > :: iterator x = state_vec_vec.begin(); x != state_vec_vec.end(); x++ ) {
			for ( vector<bool> :: iterator y = (*x).begin(); y != (*x).end(); y++ )
				sstream << *y;
			sstream << endl;
		}
		sstream << "stream" << endl;
		for ( vector< vector<bool> >::iterator x = stream_vec_vec.begin(); x != stream_vec_vec.end(); x++ ) {
			for ( vector<bool>::iterator y = (*x).begin(); y != (*x).end(); y++ )
				sstream << *y;
			sstream << endl;
		}
		known_sat_sample_file.close(); known_sat_sample_file.clear();
		known_sat_sample_file.open(known_sat_sample_file_name, ios_base::out);
		known_sat_sample_file << sstream.rdbuf();
		delete S;
	}
	else {
		string str;
		getline(known_sat_sample_file,str);
		cout << "reading state and stream from file" << endl;
		bool isState = false, isStream = false;
		do {
			if( str == "state" ) {
				cout << "state string found" << endl;
				isState = true;
			}
			else if ( str == "stream" ) {
				cout << "stream string found" << endl;
				isState = false;
				isStream = true;
			}
			else {
				if ( isState ) {
					for ( unsigned i=0; i < str.size(); i++ )
						state_vec.push_back( str[i] == '1' ? true : false );
					state_vec_vec.push_back( state_vec );
					state_vec.clear();
				}
				else if ( isStream ) {
					for ( unsigned i=0; i < str.size(); i++ )
						stream_vec.push_back( str[i] == '1' ? true : false );
					stream_vec_vec.push_back( stream_vec );
					stream_vec.clear();
				}
			}
		} while( getline(known_sat_sample_file, str ) );
		cout << "state_vec_vec.size() "  << state_vec_vec.size()  << endl;
		cout << "stream_vec_vec.size() " << stream_vec_vec.size() << endl;
	}
	cout << endl;
	known_sat_sample_file.close();
	
	if (isPlainText) // return size of input vectors without plain text data
		for (auto &x : state_vec_vec)
			x.resize(core_len);
}

string MPI_Base::MakeSolverLaunchString( string solver_name, string cnf_name, double maxtime_solving_time )
{
	string time_limit_str, result_str;
	stringstream sstream;
	sstream << max_solving_time;
	string maxtime_seconds_str = sstream.str();
	sstream.clear(); sstream.str("");
	// minisat and solvers with same time limit paremeter
	if ( solver_name.find( "cryptominisat" ) != string::npos )
		time_limit_str = "";
	else if ( ( solver_name.find( "minisat" ) != string::npos ) || 
		 ( solver_name.find( "Minisat" ) != string::npos ) || 
		 ( solver_name.find( "MiniSat" ) != string::npos ) ||
		 ( solver_name.find( "minigolf" ) != string::npos ) ||
		 ( solver_name.find( "mipisat" ) != string::npos ) ||
		 ( solver_name.find( "minitsat" ) != string::npos ) ||
		 ( solver_name.find( "rokk" ) != string::npos ) ||
		 ( solver_name.find( "sinn" ) != string::npos ) ||
		 ( solver_name.find( "zenn" ) != string::npos ) ||
		 ( solver_name.find( "SWDiA5BY" ) != string::npos ) ||
		 ( solver_name.find( "ClauseSplit" ) != string::npos ) 
		 )
	{
		time_limit_str =  "-cpu-lim=";
	}
	else if ( ( solver_name.find( "lingeling" ) != string::npos ) && 
		      ( solver_name.find( "plingeling" ) == string::npos ) ) {
		time_limit_str = "-t ";
	}
	else {
		time_limit_str = "";
		//cerr << "Unknown solver in system calling mode: " << solver_name << endl;
		//MPI_Abort( MPI_COMM_WORLD, 0 );
	}
	// glucose can't stop in time
	/*if ( solver_name.find( "glucose" ) != string::npos ) {
		cout << "glucose detected" << endl;
		result_str = "-cpu-lim=";
	}*/
	/*else if ( solver_name.find( "plingeling" ) != string::npos ) {
		//cout << "pingeling detected" << endl;
		result_str = "-nof_threads ";
		result_str += nof_threads_str;
		result_str += " -t ";
	}
	else if ( solver_name.find( "trengeling" ) != string::npos ) {
		//cout << "treengeling detected" << endl;
		//result_str = "-t " + "11" + nof_threads_str;

	if ( solver_name.find( "dimetheus" ) != string::npos )
		result_str += " -formula";
	}
	
	if ( time_limit_str == "" ) {
		cout << "unknown solver detected. using timelimit" << endl;
		result_str = "./timelimit -t " + maxtime_seconds_str + " -T 1 " + "./" + solvers_dir + "/" + solver_name;
	}*/
	
	if ( time_limit_str != "" )
		result_str = solver_name + " " + time_limit_str + maxtime_seconds_str + " " + cnf_name;
	else
		result_str = solver_name + " " + cnf_name;
	
	return result_str;
}
