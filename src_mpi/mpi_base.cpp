#include "mpi_base.h"

#pragma warning( disable : 4996 )

const int    MAX_WORD_LENGTH			   = 64;
const int    MAX_LINE_LENGTH               = 524288;
const int    MAX_LINE_LENGTH_2             = 8192;
const int    MEDIUM_STRING_LEN             = 4096;
const double TRANSP_COAST                  = 0.000001;
const int    NUM_KEY_BITS                  = 64;
//const int MAX_VAR_IN_CNF		        = 65536;
//const int MAX_SIZE_OF_FILE	        = 42949672;
//const int MAX_CNF_IN_FOLDER		    = 65536;
//const int MAX_STEP_CNF_COUNT	        = 30;
//const int LITTLE_STRING_LEN           = 16;

boost::random::mt19937 gen;

//=================================================================================================
// Constructor/Destructor:

MPI_Base :: MPI_Base( ) :
	sat_count            ( 0 ),
	corecount			 ( 0 ),
	solver_type          ( 4 ),
	core_len             ( 0 ),    
	koef_val             ( 16 ),
	schema_type          ( "0" ),
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
	IsHardProblem        ( 0 ),
	sort_type            ( 0 ),
	known_point_file_name ( "known_point" ),
	known_assumptions_file_name ( "known_assumptions" ),
	IsSolveAll           ( false ),
	IsPredict            ( false ),
	max_solving_time     ( 0 ),
	max_nof_restarts     ( 0 ),
	keybit_count         ( 4 ),
	rslos_table_name     ( "" ),
	IsFileAssumptions    ( false ),
	assumptions_string_count ( 0 ),
	activity_vec_len	 ( 0 )
{
	for ( unsigned i = 0; i < FULL_MASK_LEN; i++ )
		full_mask[i] = part_mask[i] = 0;
	gen.seed( static_cast<unsigned>(std::time(0)) );
}

MPI_Base :: ~MPI_Base( )
{
	/*int i = 0;
	if ( IsSAT ) 
		delete[] b_SAT_set_array; // allocated in ReadIntCNF
	
	if ( var_count > 0 )
	{
		for( i = 0; i < var_count; i++ )
			delete[] lits_clause_array[i];

		delete[] lits_clause_array; // allocated in ReadIntCNF
		delete[] lits_clause_lengths; // allocated in ReadIntCNF
	}

	if ( clause_count > 0 )
	{
		for( i = 0; i < clause_count; i++ )
			delete[] clause_array[i];

		delete[] clause_array; // allocated in ReadIntCNF
		delete[] clause_lengths; // allocated in ReadIntCNF
	}*/
}

static inline double cpuTime( void ) 
{
    return ( double )clock( ) / CLOCKS_PER_SEC; 
}

boost::dynamic_bitset<> MPI_Base :: IntVecToBitset( unsigned bitset_len, vector<int> &int_vec )
{
	boost::dynamic_bitset<> bs( bitset_len );
	for ( unsigned i=0; i<int_vec.size(); i++ )
		bs.set( int_vec[i]-1 ); // set element to 1
	return bs;
}

vector<int> MPI_Base :: BitsetToIntVec( boost::dynamic_bitset<> &bs )
{
	vector<int> vec_int;
	for ( unsigned i=0; i<bs.size(); i++ )
		if ( (int)bs[i] == 1 )
			vec_int.push_back( (int)(i+1) );
	return vec_int;
}

//---------------------------------------------------------
void MPI_Base :: shl64( unsigned long long int &val_for_left_shift, unsigned int bit_count )
{
	unsigned int val1, val2; 
	if ( ( bit_count > 30 ) && ( bit_count < 61 ) ) {
		val1 = 30;
		val2 = bit_count - val1;
		val_for_left_shift =  ( unsigned long long int )( 1 << val1 );
		val_for_left_shift *= ( unsigned long long int )( 1 << val2 );
	}
	else if ( bit_count < 31 )
		val_for_left_shift =  ( unsigned long long int )( 1 << bit_count );
	else
		cout << "\n bit_count " <<  bit_count << " is too large ";
}

//---------------------------------------------------------
void MPI_Base :: equalize_arr( unsigned int arr1[FULL_MASK_LEN], unsigned int arr2[FULL_MASK_LEN] )
{
	// Make arr1 be equal to arr2
	for ( int i = 0; i < FULL_MASK_LEN; i++ )
		arr1[i] = arr2[i];
}

void MPI_Base :: PrintVector( vector<int> &vec )
{
	for ( unsigned i=0; i < vec.size(); i++ )
		cout << vec[i] << " ";
	cout << endl;
}

int MPI_Base :: getdir( string dir, vector<string> &files )
{
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << endl << "Error in opening " << dir;
        return 1;
    }
    while ((dirp = readdir(dp)) != NULL) 
	{ files.push_back(string(dirp->d_name)); }
    closedir(dp);
    return 0;
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
// make array of differenr pseudorandom values
	if ( verbosity > 0 )
		cout << "MakeUniqueRandArr() started" << endl;

	if ( max_rand_val < rand_arr_len )
		max_rand_val = rand_arr_len;
	unsigned long long rand_numb;
	rand_arr.resize( rand_arr_len );
	
	bool IsOldValue;
	for ( unsigned i = 0; i < rand_arr_len; i++ ) {
		do { // if value is not unique get value again
			rand_numb = uint_rand();
			rand_numb %= max_rand_val;
			IsOldValue = false;
			for ( unsigned k = 0; k < i; k++ ) {
				if ( rand_numb == rand_arr[k] ) {
					IsOldValue = true;
					break;
				}
			}
		} while ( IsOldValue );
		rand_arr[i] = rand_numb; // new values
	}
}

// Make full_mask and part_mask for sending from order that is set by var choose array
//---------------------------------------------------------
bool MPI_Base :: GetMainMasksFromVarChoose( vector<int> &var_choose_order )
{
	if ( verbosity > 1 ) {
		cout << "GetMainMasksFromVarChoose() start with var_choose_order " << endl;
		for ( unsigned i = 0; i < var_choose_order.size(); i++ )
			cout << var_choose_order[i] << " ";
		cout << endl;
	}

	for ( unsigned i = 0; i < FULL_MASK_LEN; i++ ) {	
		full_mask[i] = 0;
		part_mask[i] = 0;
	}

	unsigned cur_uint_ind, var_index;
	for( unsigned i = 0; i < full_mask_var_count; i++ ) {
		var_index = var_choose_order[i] - 1;
		cur_uint_ind = var_index / UINT_LEN;
		full_mask[cur_uint_ind + 1] += 1 << ( var_index % UINT_LEN );
	}

	// get first full_mask_var_count vars from  array var_choose_order
	for( unsigned i = 0; i < part_mask_var_count; i++ ) {
		var_index = var_choose_order[i] - 1;
		cur_uint_ind = var_index / UINT_LEN;
		part_mask[cur_uint_ind + 1] += 1 << ( var_index % UINT_LEN );
	}

	// full_mask[0] is mask of existing all 32 uint values
	for( unsigned i = 1; i < FULL_MASK_LEN; i++ ) { // fill full_mask[0]
		if ( full_mask[i] )
			full_mask[0] += 1 << ( i-1 );
	}

	for( unsigned i = 1; i < FULL_MASK_LEN; i++ ) { // fill part_mask[0]
		if ( part_mask[i] )
			part_mask[0] += 1 << ( i-1 );
	}

	if ( ( !full_mask[0] ) && ( full_mask_var_count != 0 ) ) {
		cout << "Error. full_mask[0] == 0. full_mask_var_count != 0";
		return false;
	}

	if ( verbosity > 1 )
		cout << "GetMainMasksFromVarChoose() end" << endl;

	return true;
}

unsigned new_count_ones( unsigned num )
{
	unsigned cnt = 0;
	unsigned i;
	for ( i = 0; i < UINT_LEN; i++ ) {
		if ( num & 1 ) 
			cnt++;
		num >>= 1;
	}
	return cnt;
}

void MPI_Base :: MakeAssignsFromFile( int current_task_index, vec< vec<Lit> > &dummy_vec )
{
	ifstream in;
	string str;
	//in.open( rslos_table_name.c_str() );
	in.open( known_assumptions_file_name.c_str() );
	if ( !in.is_open() ) {
		cerr << "Error. !in.is_open(). file name " << known_assumptions_file_name << endl;
		exit;
	}

	int cur_var_ind, intval, k = 0;
	// int rslos_num = 1; 
	unsigned strings_passed = 0;
	unsigned basic_batch_size = floor( (double)assumptions_string_count / (double)all_tasks_count );
	// calculate count of bathes with additional size (+1)
	unsigned batch_addit_size_count = assumptions_string_count - basic_batch_size*all_tasks_count;
	unsigned cur_batch_size = basic_batch_size;
	if ( current_task_index < batch_addit_size_count )
		cur_batch_size++;
	// skip unuseful strings
	unsigned previous_tasks_count = current_task_index*basic_batch_size;
	if ( current_task_index < batch_addit_size_count )
		previous_tasks_count += current_task_index; // add some 1 to sum
	else
		previous_tasks_count += batch_addit_size_count;
	while ( strings_passed < previous_tasks_count ) {
		getline( in, str );
		strings_passed++;
	}
	cout << "current_task_index "       << current_task_index << endl;
	cout << "all_tasks_count "          << all_tasks_count << endl;
	cout << "previous_tasks_count "     << previous_tasks_count << endl;
	cout << "assumptions_string_count " << assumptions_string_count << endl;
	cout << "basic_batch_size "         << basic_batch_size  << endl;
	cout << "cur_batch_size "           << cur_batch_size << endl;
	cout << "strings_passed "           << strings_passed  << endl;
	
	dummy_vec.growTo( cur_batch_size );
	// reading values from file
	for ( unsigned i=0; i < cur_batch_size; i++ ) {
		if ( !getline( in, str ) ) {
			cerr << "Error. !getline( in, str )" << endl;
			exit;
		}
		//rslos_num = 1;
		//k = 0;
		if ( str.size() < var_choose_order.size() ) {
			cerr << "Error. str.size() < var_choose_order.size()" << endl;
			exit;
		}
		for ( unsigned j=0; j < var_choose_order.size(); j++ ) {
			cur_var_ind = var_choose_order[j] - 1;
			intval = j; //intval = keybit_count*rslos_num - k - 1;
			if ( str[intval] == '1' ) {
				dummy_vec[i].push( mkLit( cur_var_ind ) );
				//cout << cur_var_ind << " ";
			} else {
				dummy_vec[i].push( ~mkLit( cur_var_ind ) );
				//cout << "-" << cur_var_ind << " ";
			}
			/*if ( k == keybit_count - 1 ) {
				rslos_num++;
				k = 0;
			} else k++;*/
		}
	}
	in.close();
}

void MPI_Base :: MakeAssignsFromMasks( unsigned full_mask[FULL_MASK_LEN], 
									   unsigned part_mask[FULL_MASK_LEN], 
									   unsigned value[FULL_MASK_LEN],
									   vec< vec<Lit> > &dummy_vec )
{
// for predict with minisat2.2. convert masks to vector of Literals
	unsigned mask, range = 0, range_mask_ind;
	int cur_var_ind;
	unsigned long long range_val_count,
			           range_mask, 
			           lint;
	bool IsPositiveLiteral, IsAddingLiteral;
	Lit new_lit;

	for ( unsigned i = 1; i < FULL_MASK_LEN; i++ )
		range += new_count_ones( full_mask[i] ^ part_mask[i] );

	shl64( range_val_count, range );
	dummy_vec.growTo( range_val_count );

	for ( lint = 0; lint < range_val_count; lint++ ) {
		range_mask = 1;
		range_mask_ind = 1;

		for ( unsigned i = 1; i < FULL_MASK_LEN; i++ ) {
			for ( unsigned j = 0; j < UINT_LEN; j++ ) {
				mask = ( 1 << j );
				cur_var_ind = ( i-1 ) * UINT_LEN + j;
				IsAddingLiteral = false;
				if ( part_mask[i] & mask ) {
					IsPositiveLiteral = ( value[i] & mask ) ? true : false;
					IsAddingLiteral = true;
				}
				else if ( full_mask[i] & mask ) {
					IsPositiveLiteral = ( lint & range_mask ) ? true : false;
 					shl64( range_mask, range_mask_ind );
					range_mask_ind++;
					IsAddingLiteral = true;
				}
				if ( !IsAddingLiteral ) continue;

				new_lit = IsPositiveLiteral ? mkLit( cur_var_ind ) : ~mkLit( cur_var_ind );
				dummy_vec[lint].push( new_lit );
			}
		}
	}
}

// Get values for sending using order by var choose array
//---------------------------------------------------------
bool MPI_Base :: GetValuesFromVarChoose( unsigned long long int &part_var_power,
										 unsigned int **&values_arr )
{
	unsigned int cur_uint_ind = 0;
	int i, j;
	
	unsigned long long int lint;
	unsigned int value_index;
	unsigned int  mask;
	unsigned long long int mask2;
	for( lint = 0; lint < part_var_power; lint++ ) {
		value_index = 0;
		for ( i = 1; i < FULL_MASK_LEN; i++ ) {
			for( j = 0; j < UINT_LEN; j++ ) {
				mask  = ( 1 << j );
				shl64( mask2, value_index );
				if ( part_mask[i] & mask ) { // if part_mask bit is 1	
					if ( lint & mask2 ) // if bit is 1
						values_arr[lint][i] += mask;
					value_index++;
				}
			} // for( j = 0; j < UINT_LEN; j++ )
		}
	} // for( lint = 0; lint < part_var_power; lint++ )

	// set values_arr[lint][0] and values_arr[lint][i]
	for( lint = 0; lint < part_var_power; lint++ ) {
		for( i = 1; i < FULL_MASK_LEN; i++ ) { // fill part_mask[0]
			if ( values_arr[lint][i] )
				values_arr[lint][0] += 1 << ( i-1 );
		}
	}

	return true;
}

//---------------------------------------------------------
bool MPI_Base :: MakeStandartMasks( unsigned long long int &part_var_power, 
									unsigned int **&values_arr )
{		
	if ( !GetMainMasksFromVarChoose( var_choose_order ) ) { 
		cout << "\n Error in GetMainMasksFromVarChoose" << endl; 
		return false; 
	}
	cout << "\n Correct end of GetMasksFromVarChoose" << endl;
		
	if ( !GetValuesFromVarChoose( part_var_power, values_arr ) ) { 
		cout << "\n Error in GetValuesFromVarChoose" << endl; 
		return false; 
	}
	cout << "\n Correct end of GetValuesFromVarChoose" << endl;
	
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

	ifstream known_assumptions_file( known_assumptions_file_name.c_str() );
	if ( known_assumptions_file.is_open() ) {
		IsFileAssumptions = true;
		string str;
		while ( getline( known_assumptions_file, str ) ) {
			if ( str.size() > 1 ) // skip empty strings
			assumptions_string_count++;
		}
		known_assumptions_file.close();
	}
	
	// if file with decomposition set exists
	ifstream known_point_file( known_point_file_name.c_str() );
	if ( known_point_file.is_open() ) {
		getline( known_point_file, str );
		sstream << str;
		while ( sstream >> val )
			var_choose_order.push_back( val );
		full_mask_var_count = var_choose_order.size();
		sstream.str( "" ); sstream.clear( );
		known_point_file.close();
		schema_type = "known_point";
	}
	
	if ( ( schema_type == "rslos_end" ) && ( rslos_lengths.size() > 0 ) ) {
		IsFileAssumptions = true;
		known_assumptions_file_name = rslos_table_name;
		cout << "schema_type " << schema_type << endl;
		int keybit_table = 4; // TODO change to variable
		int prev_lengths_sum = 0;
		full_mask_var_count = rslos_lengths.size() * keybit_table; 
		for ( unsigned i=0; i < rslos_lengths.size(); ++i ) {
			for ( int j=0; j < keybit_table; ++j )
				var_choose_order.push_back( prev_lengths_sum + rslos_lengths[i] - j );
			prev_lengths_sum += rslos_lengths[i];
		}
	}

	cout << "full_mask_var_count " << full_mask_var_count << endl;
	var_choose_order.resize( full_mask_var_count );

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
	else if ( schema_type == "2" ) {
			// Literal count
			// init before sorting
			for ( unsigned  i = 0; i < core_len; i++ )
				var_choose_order[i] = i + 1;
			var_literal_count_weights = new unsigned[core_len];
			for ( unsigned i = 0; i < core_len; i++ )
				var_literal_count_weights[i] = lits_clause_lengths[i*2] + lits_clause_lengths[i*2 + 1];
			/*for ( i = 0; i < core_len; i++ )
				cout << var_literal_count_weights[i] << endl;
			cout << endl;*/
			/*for ( int i = core_len - 1; i > -1 ; --i) // bubble sort
				for ( int j = 0; j < i; j++ )
					if ( var_literal_count_weights[j] < var_literal_count_weights[j + 1] ) {
						tmp = var_literal_count_weights[j];
						var_literal_count_weights[j] = var_literal_count_weights[j + 1];
						var_literal_count_weights[j + 1] = tmp;
						tmp = var_choose_order[j];
						var_choose_order[j] = var_choose_order[j + 1];
						var_choose_order[j + 1] = tmp;
					}*/
			/*cout << "var_literal_count_weights" << endl;
			for ( i = 0; i < core_len; i++ )
				cout << var_literal_count_weights[i] << " ";
			cout << endl;*/
			delete[] var_literal_count_weights;
		}
	else if ( schema_type == "3" ) { // Jeroslaw-Wang heruistic, sum(1/2^len(clauses))
			// init before sorting
			for ( unsigned i = 0; i < core_len; i++ )
				var_choose_order[i] = i + 1;
			var_jeroslaw_count_weights = new double[core_len];

			for ( unsigned i = 0; i < core_len; i++ ) {
				var_jeroslaw_count_weights[i] = 0;
				for ( unsigned j = 0; j < lits_clause_lengths[i*2]; j++ ) // sum for positiv literal
					var_jeroslaw_count_weights[i] += 1/( pow( ( double)2, clause_lengths[lits_clause_array[i*2][j]] ) );
				for ( unsigned j = 0; j < lits_clause_lengths[i*2 + 1]; j++ ) // // sum for negative literal
					var_jeroslaw_count_weights[i] += 1/( pow( ( double)2, clause_lengths[lits_clause_array[i*2 + 1][j]] ) );
			}

			/*for ( int i = core_len - 1; i > -1 ; --i) // bubble sort
				for ( int j = 0; j < i; j++)
					if ( var_jeroslaw_count_weights[j] < var_jeroslaw_count_weights[j + 1] ) {
						double_tmp = var_jeroslaw_count_weights[j];
						var_jeroslaw_count_weights[j] = var_jeroslaw_count_weights[j + 1];
						var_jeroslaw_count_weights[j + 1] = double_tmp;
						tmp = var_choose_order[j];
						var_choose_order[j] = var_choose_order[j + 1];
						var_choose_order[j + 1] = tmp;
					}*/

			delete[] var_jeroslaw_count_weights;
	}
	else if ( schema_type == "4" ) {
	// Implicant count heruistic, for z =  count of literals that are not in same clause
			// init before sorting
			for ( unsigned i = 0; i < core_len; i++ )
				var_choose_order[i] = i + 1;
			var_implicant_count_weights = new unsigned[core_len];
			
			for ( unsigned i = 0; i < core_len; i++ ) { // sum for positiv literal
				var_implicant_count_weights[i] = lit_count - 1; // minus contrar literal
				for ( unsigned j = 0; j < lits_clause_lengths[i*2]; j++ ) // for all clauses which has literal # (i + 1)*2
					var_implicant_count_weights[i] -= clause_lengths[lits_clause_array[i*2][j]];
			}

			for ( unsigned i = 0; i < core_len; i++ ) // sum for negativ literal
			{
				var_implicant_count_weights[i] += lit_count - 1; // minus contrar literal
				for ( unsigned j = 0; j < lits_clause_lengths[i*2 + 1]; j++ ) // for all clauses which has literal # (i + 1)*2 + 1
					var_implicant_count_weights[i] -= clause_lengths[lits_clause_array[i*2 + 1][j]];
			}

			/*cout << "var_implicant_count_weights_weights" << endl;
			for ( i = 0; i < core_len; i++ )
				cout << var_implicant_count_weights[i] << " ";
			cout << endl;*/

			/*for ( int i = core_len - 1; i > -1 ; --i) // bubble sort
				for ( int j = 0; j < i; j++)
					if ( var_implicant_count_weights[j] < var_implicant_count_weights[j + 1] ) {
						tmp = var_implicant_count_weights[j];
						var_implicant_count_weights[j] = var_implicant_count_weights[j + 1];
						var_implicant_count_weights[j + 1] = tmp;
						tmp = var_choose_order[j];
						var_choose_order[j] = var_choose_order[j + 1];
						var_choose_order[j + 1] = tmp;
					}*/
			delete[] var_implicant_count_weights;
	}

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

	clause_lengths      = new int[clause_count];
	lits_clause_lengths = new unsigned[lit_count];
	for ( i = 0; i < lit_count; i++ )
		lits_clause_lengths[i] = 0;

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
					// convert value of literal to positiv int
					lit_positiv_val = abs( atoi( word_str.c_str( ) ) ); 
					if ( lit_positiv_val > current_var_count )
						current_var_count = lit_positiv_val;

					lit_val = atoi( word_str.c_str( ) ); // word -> lit value
					if ( lit_val < 0 ) {
						lit_val = -lit_val; 
						sign = 1;
					}
					else
						sign = 0;

					// add values to attay of literal's clauses relations
					val = ( lit_val << 1 ) + sign; // literal value, 1 -> 2, -1 -> 3, 2 -> 4, -2 -> 5
					lit_index = val - 2;
					lits_clause_lengths[lit_index] += 1; 
				}
			} // for ( i = 0; i < line_str.length( ) - 1; i++ )

			if ( current_lit_count )
			{
				current_clause_count++;
				clause_lengths[current_clause_count - 1] = current_lit_count;
			}
		}
	}

	/*for ( i = 0; i < 20; i++ )
		cout << clause_lengths[i] << " " << endl;
	cout << endl;
	for ( i = 0; i < 20; i++ )
		cout << lits_clause_lengths[i] << " " << endl;*/

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
	unsigned i,
		input_var_num,
		first_obj_var;
	int lit_val;
	string line_str, 
		   word_str;
	bool IncorrectLine;

	if ( !ReadVarCount( ) ) {
		cerr << "Error in ReadVarCount" << endl;
		return false;
	}
	if ( verbosity > 0 )
		cout << "Success of ReadVarCount()" << endl;

	b_SAT_set_array = new int[var_count];

	clause_array      = new int*[clause_count];
	lits_clause_array = new int*[lit_count]; // array of indexes of clauses for literals

	int *lits_clause_current = new int[lit_count];

	for ( i = 0; i < clause_count; i++ )
		clause_array[i]      = new int[clause_lengths[i]];
	for ( i = 0; i < lit_count; i++ ) {
		lits_clause_array[i] = new int[lits_clause_lengths[i]];
		lits_clause_current[i] = 0;
	}

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
			
			if ( str2 == "rslos" ) {
				while ( sstream >> intval )
					rslos_lengths.push_back( intval );
				cout << "rslos_count " << rslos_lengths.size() << endl;
				cout << "rslos lens: ";
				for ( unsigned i = 0; i < rslos_lengths.size(); i++ )
					cout << rslos_lengths[i] << " ";
				cout << endl;
				
				continue;
			}

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
				    Is_InpVar = true;
				    continue;
				}
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
								cout << "\n ***Error. String c obj vars ... contains too many values ";
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
			for ( i = 0; i < line_str.length( ) - 1; i++ ) {
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
					// fill attay of literal's clauses indexes
					int lit_index = val - 2;
					lits_clause_array[lit_index][lits_clause_current[lit_index]] = ( int )current_clause_count;
					lits_clause_current[lit_index]++;
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

	delete[] lits_clause_current;
	main_cnf.close( );

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
		str_output_cur_ind = 0, 
		k;
	unsigned i , answer_var_count;
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
	//cout << "var_count " << var_count << endl;
	lit_SAT_set_array.resize( var_count );
	for ( i = 0; i < lit_SAT_set_array.size(); i++ ) {
		lit_SAT_set_array[i] = ( ( i + 1 ) << 1 ) + (b_SAT_set_array[i] ? 0 : 1);
		//cout << lit_SAT_set_array[i];
	}
	cout << endl;
	// if SAT set exist then check it
	//cout << "Before CheckSATset()" << endl;
	if ( !CheckSATset( lit_SAT_set_array ) ) {
		cerr << "Error in checking of SAT set" << endl;
		return false;
	}
	if ( verbosity > 0 )
		cout << "CheckSATset() done" << endl;
	
	// open file for writing of answer in zchaff style
	answer_file_name = "sat_sets_";
	answer_file_name += input_cnf_name;
	
	k = 1;
	for( i = 1; i < ( int )answer_file_name.length( ); i++ ) {
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
	
	// put only core vars if such was specified
	if ( !core_len ) answer_var_count = var_count;
	else answer_var_count = core_len;
	
	sstream << "SAT" << endl;
	for ( unsigned i = 0; i < var_count; i++ )
		//sstream << (b_SAT_set_array[i] ? i + 1 : -( i + 1 )) << " ";
		sstream << b_SAT_set_array[i];
	sstream << endl; 
	answer_file << sstream.rdbuf( );
	answer_file.close( );
	lit_SAT_set_array.clear();
	
	return true;
}

//---------------------------------------------------------
bool MPI_Base :: SortValuesDecrease( unsigned int range, 
									 unsigned int *&sorted_index_array )
{
// Sorting of array {0, ... , 2^range - 1 }
// by decreasing of weight (from 1-vector to 0-vector)
// bits - array of "1" count in 0 .. 255
	unsigned int bits[] = {0,1,1,2,1,2,2,3,1,
      2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,
      2,3,3,4,3,4,4,5,1,2,2,3,2,3,3,4,2,3,
      3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,
      4,5,5,6,1,2,2,3,2,3,3,4,2,3,3,4,3,4,
      4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
      2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,
      4,5,4,5,5,6,4,5,5,6,5,6,6,7,1,2,2,3,
      2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,
      4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,
      3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,
      5,6,5,6,6,7,2,3,3,4,3,4,4,5,3,4,4,5,
      4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,
      6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
      4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8};

	// part_var_power - length of binary sequences for sorting
	unsigned int j = 0,
				 bit_weight = 0,
				 k,
				 val2;
	int i, 
		table_val_count = 0, 
		bit_value_count = 0;
	unsigned int *index_array_lengths,
				 *sort_index_arr, 
				 **index_array;

	table_val_count = ( 1 << range );
	bit_value_count = range + 1;
	
	index_array_lengths = ( unsigned int* )malloc( bit_value_count * sizeof( unsigned int ) );
	sort_index_arr = ( unsigned int* )malloc( bit_value_count * sizeof( unsigned int ) );
	index_array = ( unsigned int** )malloc( bit_value_count * sizeof( unsigned int* ) );

	for ( i = 0; i < bit_value_count; i++ )
		index_array_lengths[i] = 0;
	
	// at first compute lengths of arrays, they are same for any kind of input values
	for ( i = 0; i < table_val_count; i++ )
	{
		// get count of "1" in 32-bit value by shifting 3 times
		bit_weight = bits[( unsigned char )( i )] + bits[( unsigned char )( i >> 8 )] +
		             bits[( unsigned char )( i >> 16 )] + bits[( unsigned char )( i >> 24 )];
		index_array_lengths[bit_weight] += 1;
		//index_array[bit_weight][index_array_lengths[bit_weight] - 1] = i;
	}
	
	// allocate array of indexes
	for ( i = 0; i < bit_value_count; i++ )
	{
		if ( index_array_lengths[i] )
			index_array[i] = ( unsigned int* )malloc( index_array_lengths[i] * sizeof( unsigned int ) );
		index_array_lengths[i] = 0; // null array of lengths
	}
	
	// fill array of indexes
	for ( i = 0; i < table_val_count; i++ )
	{
		// by shifting 3 times get count of "1" in 32-bit value
		bit_weight = bits[( unsigned char )( i )] + bits[( unsigned char )( i >> 8 )] +
		             bits[( unsigned char )( i >> 16 )] + bits[( unsigned char )( i >> 24 )];
		// again fill array of lengths
		index_array_lengths[bit_weight] += 1;
		index_array[bit_weight][index_array_lengths[bit_weight] - 1] = i;
	}

	k = 0;
	// fill array with sorted values
	for ( int i = bit_value_count - 1; i > -1; --i )
	{
		val2 = index_array_lengths[i];
		for ( j = 0; j < val2; j++ ) {
			sorted_index_array[k] = index_array[i][j];
			k++;
		}
	}
	
	// deallocate arrays
	free( sort_index_arr );

	for ( i = 0; i < bit_value_count; i++ ) {
		if ( index_array_lengths[i] )
			free( index_array[i] );
	}
	
	free( index_array_lengths );
	free( index_array );
	
	return true;
}

unsigned MPI_Base ::  uint_rand() {
	boost::random::uniform_int_distribution<uint32_t> dist;
	return dist(gen);
}

void MPI_Base :: AddSolvingTimeToArray( ProblemStates cur_problem_state, double cnf_time_from_node, double *solving_times )
{
	// solving_times[0]  == min
	// solving_times[1]  == max
	// solving_times[2]  == med
	// solving_times[3]  == sat
	switch( cur_problem_state ){ 
		case Solved :
			if ( cnf_time_from_node < solving_times[0] )
				solving_times[0] = cnf_time_from_node;
			if ( cnf_time_from_node > solving_times[1] ) {
				solving_times[1] = cnf_time_from_node;
			}
			else if ( cnf_time_from_node < 0.0001 ) solving_times[5]++;
			else if ( cnf_time_from_node < 0.001  ) solving_times[6]++;
			else if ( cnf_time_from_node < 0.01 )   solving_times[7]++;
			else if ( cnf_time_from_node < 0.1  )   solving_times[8]++;
			else if ( cnf_time_from_node < 1    )   solving_times[9]++;
			else if ( cnf_time_from_node < 10   )   solving_times[10]++;
			else if ( cnf_time_from_node < 100  )   solving_times[11]++;
			else if ( cnf_time_from_node < 1000 )   solving_times[12]++;
			else solving_times[13]++;
			break;
		case SolvedOnPreprocessing : 
			solving_times[4]++;
			break;
		case Interrupted : 
			solving_times[14]++;
			break;
	}
}

bool MPI_Base :: SolverRun( Solver *&S, unsigned int *full_mask, unsigned int *part_mask, 
							unsigned int *value, int &process_sat_count, int &current_obj_val,
							double &cnf_time_from_node, double *solving_times, 
							int current_task_index )
{
// Run needed solver
	process_sat_count = 0;
	
	ofstream file_class_prep, file_class1, file_class2, file_class3, 
		     file_class4, file_class5, file_class_sat;
	stringstream sstream;
	vec< vec<Lit> > dummy_vec;
	lbool ret;
	double total_time = 0;
	unsigned current_tasks_solved = 0;
	int result;
	ProblemStates cur_problem_state;
	
	if ( verbosity > 1 )
		cout << "start SolverRun()" << endl;

	solving_times[0] = 1 << 30; // start min len
	for ( unsigned i = 1; i < SOLVING_TIME_LEN; i++ )
		solving_times[i] = 0;
	
	if ( solver_type == 1 ) {
		//printf( "\n Start of dminisat_solve" );
		// param 0 means flag IsPredict == 0
		if ( !dminisat_solve( input_cnf_name, full_mask, part_mask, value,
							  solver_type, core_len, start_activity, 
							  &process_sat_count, &b_SAT_set_array, sort_type,
							  &cnf_time_from_node, 0, IsHardProblem ) )
		{ printf( "\n Error in dminisat_solve" ); }
	} 
	else if ( solver_type == 4 ) {	
		if ( IsFileAssumptions ) // if assumptions in file 
			MakeAssignsFromFile( current_task_index, dummy_vec );
		else
			MakeAssignsFromMasks( full_mask, part_mask, value, dummy_vec );
		
		if ( verbosity > 1 ) {
			cout << "dummy_vec size" << dummy_vec.size() << endl;
			for ( int i = 0; i < dummy_vec.size(); i++ ) {
				for ( int j=0; j < dummy_vec[i].size(); j++ )
					cout << dummy_vec[i][j].x << " ";
				cout << endl;
			}
		}

		uint64_t prev_starts, prev_conflicts, prev_decisions;
		for ( int i=0; i < dummy_vec.size(); i++ ) {
#ifndef _DEBUG
			cnf_time_from_node = MPI_Wtime( );
#endif
			// save current state to check differences
			prev_starts    = S->starts;
			prev_conflicts = S->conflicts;
			prev_decisions = S->decisions;
			
			S->last_time = Minisat :: cpuTime();
			ret = S->solveLimited( dummy_vec[i] );
			
			//ret = S->solveLimited( dummy_vec[i], true, false ); // for SimpSolver
#ifndef _DEBUG
			cnf_time_from_node = MPI_Wtime( ) - cnf_time_from_node;
#endif
			total_time += cnf_time_from_node;

			if ( ret == l_Undef )
				cur_problem_state = Interrupted; // interrupted cause of restarts or time limit
			else if ( ( S->starts - prev_starts <= 1 ) && ( S->conflicts == prev_conflicts ) && ( S->decisions == prev_decisions ) )
				cur_problem_state = SolvedOnPreprocessing;  // solved by BCP
			else
				cur_problem_state = Solved; // just solved
			
			AddSolvingTimeToArray( cur_problem_state, cnf_time_from_node, solving_times );
			result = (ret == l_True) ? 1 : 0;
			
			if ( result ) {
				process_sat_count += result;
				cout << "process_sat_count " << process_sat_count << endl;
				if ( !solving_times[3] ) {
					solving_times[3] = cnf_time_from_node; // time of 1st SAT, write only once
					cout << " SAT time " << solving_times[3] << endl;
				}
				for ( int i=0; i < S->model.size(); i++ )
					b_SAT_set_array[i] = ( S->model[i] == l_True) ? 1 : 0 ;
						// check res file for SAT set existing
				if ( !AnalyzeSATset( ) ) {
					// is't needed to deallocate memory - MPI_Abort will do it	
					cout << "\n Error in Analyzer" << endl;
					MPI_Abort( MPI_COMM_WORLD, 0 );
					return false;
				}
				if ( !IsSolveAll )
					break;
			}
			//S->loadState();
			//delete S;
			//S->clearDB(); // we don't need to clear DB, because incremental solving is faster
		}
		
		solving_times[2] = total_time / dummy_vec.size(); // med time for current batch
		//cout << "median time in batch " << solving_times[2] << endl;

		for ( int i = 0; i < dummy_vec.size(); i++ )
			dummy_vec[i].clear();
		dummy_vec.clear();
	}
	else 
		{ cout << "\n solver_type has unknown format" << endl; return false; }

	return true;
}
