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

//=================================================================================================
// Constructor/Destructor:

MPI_Base :: MPI_Base( ) :
	IsSAT                ( false ),
	corecount			 ( 0 ),
	solver_type          ( 1 ),
	core_len             ( 0 ),    
	koef_val             ( 4 ),
	schema_type          ( 0 ),
	var_count            ( 0 ),
	clause_count         ( 0 ),
	full_mask_var_count  ( 0 ),
	corevars_activ_type  ( 1 ),
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
	verbosity			 ( 1 ),
	check_every_conflict ( 2000 )
{
	unsigned int i;
	for ( i = 0; i < FULL_MASK_LEN; i++ )
	{	
		full_mask[i] = 0;
		part_mask[i] = 0;
	}
	for ( i = 0; i < MAX_CORE_LEN; i++ )
		var_choose_order[i] = -1;
}

MPI_Base :: ~MPI_Base( )
{
	int i = 0;
	/*if ( IsSAT ) 
		delete[] b_SAT_set_array; // allocated in ReadIntCNF
	*/
	
	if ( var_count > 0 )
	{
		delete[] b_SAT_set_array;
		for( i = 0; i < var_count * 2; i++ )
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
	}
}

static inline double cpuTime( void ) 
{
    return ( double )clock( ) / CLOCKS_PER_SEC; 
}

//---------------------------------------------------------
void MPI_Base :: shl64( unsigned long long int &val_for_left_shift, unsigned int bit_count )
{
	unsigned int val1, val2; 
	if ( ( bit_count > 30 ) && ( bit_count < 61 ) )
	{
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

// Make full_mask and part_mask for sending from order that is set by var choose array
//---------------------------------------------------------
bool MPI_Base :: GetMainMasksFromVarChoose( )
{
	int i, cur_uint_ind;

	for ( i = 0; i < FULL_MASK_LEN; i++ )
	{	
		full_mask[i] = 0;
		part_mask[i] = 0;
	}

	for( i = 0; i < full_mask_var_count; i++ )
	{
		cur_uint_ind = var_choose_order[i] / UINT_LEN;
		full_mask[cur_uint_ind + 2] += 1 << ( var_choose_order[i] % UINT_LEN);
	}

	/*for( i = 0; i < part_mask_var_count; i++ )
		printf( "\n var_choose_order[i] %d", var_choose_order[i] );*/

	// get first full_mask_var_count vars from  array var_choose_order
	for( i = 0; i < part_mask_var_count; i++ )
	{
		cur_uint_ind = var_choose_order[i] / UINT_LEN;
		part_mask[cur_uint_ind + 2] += 1 << ( var_choose_order[i] % UINT_LEN );
	}

	bool notnull_uint = true;
	// full_mask[0] is mask of existing all 32 uint values
	// part_mask[1] is count (0..32) of first not null uint in 2..34 range
	for( i = 2; i < FULL_MASK_LEN; i++ ) // fill full_mask[0] and full_mask[1] 
	{
		if ( full_mask[i] ) 
		{
			full_mask[0] += 1 << ( i - 2 );
			if ( notnull_uint )
				full_mask[1]++;
		}
		else
			notnull_uint = false; 
	}

	/*
	if ( full_mask[1] > 32 )
		printf("\n Error");
	*/

	notnull_uint = true;
	for( i = 2; i < FULL_MASK_LEN; i++ ) // fill part_mask[0] and part_mask[1] 
	{
		if ( part_mask[i] ) 
		{
			part_mask[0] += 1 << ( i - 2 );
			if ( notnull_uint )
				part_mask[1]++;
		}
		else 
			notnull_uint = false;
	}

	if ( ( full_mask[1] == 0 ) && ( full_mask[0] != 0 ) ) // if 1st uint == 0
		full_mask[1] = UINT_LEN;
	
	if ( ( part_mask[1] == 0 ) && ( part_mask[0] != 0 ) ) // if 1st uint == 0
		part_mask[1] = UINT_LEN;

	if ( ( !full_mask[0] ) && ( full_mask_var_count != 0 ) )
	{
		std :: cout << "Error. full_mask[0] == 0. full_mask_var_count != 0";
		return false;
	}

	return true;
}

// Get values for sending from order that is set by var choose array
//---------------------------------------------------------
bool MPI_Base :: GetValuesFromVarChoose( unsigned long long int &part_var_power,
										 unsigned int **&values_arr )
{
	unsigned int cur_uint_ind = 0;
	int i, j;
	bool notnull_uint;
	
	unsigned long long int lint;
	unsigned int value_index;
	unsigned int  mask;
	unsigned long long int mask2;
	for( lint = 0; lint < part_var_power; lint++ )
	{
		value_index = 0;
		for ( i = 2; i < ( int )( part_mask[1] + 2 ); i++ )
		{
			for( j = 0; j < UINT_LEN; j++ )
			{
				mask  = ( 1 << j );
				shl64( mask2, value_index );
				if ( part_mask[i] & mask ) // if part_mask bit is 1
				{	
					if ( lint & mask2 ) // if bit is 1
						values_arr[lint][i] += mask;
					value_index++;
				}
			} // for( j = 0; j < UINT_LEN; j++ )
		} // for ( i = 2; i < part_mask[1] + 2; i++ )
	} // for( lint = 0; lint < part_var_power; lint++ )

	// set values_arr[lint][0] and values_arr[lint][i]
	for( lint = 0; lint < part_var_power; lint++ )
	{
		notnull_uint = true;
		for( i = 2; i < FULL_MASK_LEN; i++ ) // fill part_mask[0] and part_mask[1] 
		{
			if ( values_arr[lint][i] )
			{
				values_arr[lint][0] += 1 << ( i - 2 );
				if ( notnull_uint )
					values_arr[lint][1]++;
			}
			else
				notnull_uint = false;
		}

		if ( ( values_arr[lint][1] == 0 ) && ( values_arr[lint][0] != 0 ) )// if 1st uint == 0
			values_arr[lint][1] = UINT_LEN;
	}

	return true;
}

//---------------------------------------------------------
bool MPI_Base :: MakeStandartMasks( unsigned long long int &part_var_power, 
									unsigned int **&values_arr )
{
	if ( !MakeVarChoose( ) ) { 
		cout << "\n Error in MakeVarChoose" << endl; 
		return false; 
	}
	cout << "\n Correct end of MakeVarChoose" << endl;
		
	if ( !GetMainMasksFromVarChoose( ) ) { 
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
	int i = 0, k = 0, j = 0;
	int LFSR_1_arr_len,
		LFSR_2_arr_len,
		LFSR_3_arr_len,
		LFSR_4_arr_len,
		LFSR_5_arr_len,
		LFSR_1_mid_ind,
		LFSR_2_mid_ind,
		LFSR_3_mid_ind,
		LFSR_1_cur_ind,
		LFSR_2_cur_ind,
		LFSR_3_cur_ind,
		LFSR_4_cur_ind,
		LFSR_5_cur_ind,
		sum4_66_LFSR_1_len = 13,
		sum4_66_LFSR_2_len = 15,
		sum4_66_LFSR_3_len = 17,
		sum4_66_LFSR_4_len = 19,
		Gifford_block_len  = 8,
		tmp;
	unsigned int des_key_len = 56;
	int first_var_init = -1;
	int *var_literal_count_weights;
	int *var_implicant_count_weights;
	double double_tmp;
	double *var_jeroslaw_count_weights;
		
	switch ( schema_type )
	{
		case 0: // first vars: 0, 1, 2, ...
			for ( i = 0; i < full_mask_var_count; i++ )
				var_choose_order[i] = i;
			break;
		case 1: // missone: 0, 2, 4, ...
			for ( i = 0; i < full_mask_var_count; i++ )
				var_choose_order[i] = i*2;
			break;
		case 2: // Literal count
			// init before sorting
			for ( i = 0; i < core_len; i++ )
				var_choose_order[i] = i;
			var_literal_count_weights = new int[core_len];
			for ( i = 0; i < core_len; i++ )
				var_literal_count_weights[i] = lits_clause_lengths[i*2] + lits_clause_lengths[i*2 + 1];
			/*for ( i = 0; i < core_len; i++ )
				cout << var_literal_count_weights[i] << endl;
			cout << endl;*/
			for ( i = core_len - 1; i > -1 ; --i) // bubble sort
				for ( j = 0; j < i; j++)
					if ( var_literal_count_weights[j] < var_literal_count_weights[j + 1] )
					{
						tmp = var_literal_count_weights[j];
						var_literal_count_weights[j] = var_literal_count_weights[j + 1];
						var_literal_count_weights[j + 1] = tmp;
						tmp = var_choose_order[j];
						var_choose_order[j] = var_choose_order[j + 1];
						var_choose_order[j + 1] = tmp;
					}
			/*cout << "var_literal_count_weights" << endl;
			for ( i = 0; i < core_len; i++ )
				cout << var_literal_count_weights[i] << " ";
			cout << endl;*/
			delete[] var_literal_count_weights;
			break;
		case 3: // Jeroslaw-Wang heruistic, sum(1/2^len(clauses))
			// init before sorting
			for ( i = 0; i < core_len; i++ )
				var_choose_order[i] = i;
			var_jeroslaw_count_weights = new double[core_len];

			for ( i = 0; i < core_len; i++ )
			{
				var_jeroslaw_count_weights[i] = 0;
				for ( j = 0; j < lits_clause_lengths[i*2]; j++ ) // sum for positiv literal
					var_jeroslaw_count_weights[i] += 1/( pow( ( double)2, clause_lengths[lits_clause_array[i*2][j]] ) );
				for ( j = 0; j < lits_clause_lengths[i*2 + 1]; j++ ) // // sum for negative literal
					var_jeroslaw_count_weights[i] += 1/( pow( ( double)2, clause_lengths[lits_clause_array[i*2 + 1][j]] ) );
			}

			/*cout << "var_jeroslaw_count_weights" << endl;
			for ( i = 0; i < core_len; i++ )
				cout << var_jeroslaw_count_weights[i] << " ";
			cout << endl;*/

			for ( i = core_len - 1; i > -1 ; --i) // bubble sort
				for ( j = 0; j < i; j++)
					if ( var_jeroslaw_count_weights[j] < var_jeroslaw_count_weights[j + 1] )
					{
						double_tmp = var_jeroslaw_count_weights[j];
						var_jeroslaw_count_weights[j] = var_jeroslaw_count_weights[j + 1];
						var_jeroslaw_count_weights[j + 1] = double_tmp;
						tmp = var_choose_order[j];
						var_choose_order[j] = var_choose_order[j + 1];
						var_choose_order[j + 1] = tmp;
					}
			/*cout << "var_jeroslaw_count_weights" << endl;
			for ( i = 0; i < core_len; i++ )
				cout << var_jeroslaw_count_weights[i] << " ";
			cout << endl;*/

			delete[] var_jeroslaw_count_weights;
			break;
		case 4: // Implicant count heruistic, for z =  count of literals that are not in same clause
			// init before sorting
			for ( i = 0; i < core_len; i++ )
				var_choose_order[i] = i;
			var_implicant_count_weights = new int[core_len];
			
			for ( i = 0; i < core_len; i++ ) // sum for positiv literal
			{
				var_implicant_count_weights[i] = lit_count - 1; // minus contrar literal
				for ( j = 0; j < lits_clause_lengths[i*2]; j++ ) // for all clauses which has literal # (i + 1)*2
					var_implicant_count_weights[i] -= clause_lengths[lits_clause_array[i*2][j]];
			}

			for ( i = 0; i < core_len; i++ ) // sum for negativ literal
			{
				var_implicant_count_weights[i] += lit_count - 1; // minus contrar literal
				for ( j = 0; j < lits_clause_lengths[i*2 + 1]; j++ ) // for all clauses which has literal # (i + 1)*2 + 1
					var_implicant_count_weights[i] -= clause_lengths[lits_clause_array[i*2 + 1][j]];
			}

			/*cout << "var_implicant_count_weights_weights" << endl;
			for ( i = 0; i < core_len; i++ )
				cout << var_implicant_count_weights[i] << " ";
			cout << endl;*/

			for ( i = core_len - 1; i > -1 ; --i) // bubble sort
				for ( j = 0; j < i; j++)
					if ( var_implicant_count_weights[j] < var_implicant_count_weights[j + 1] )
					{
						tmp = var_implicant_count_weights[j];
						var_implicant_count_weights[j] = var_implicant_count_weights[j + 1];
						var_implicant_count_weights[j + 1] = tmp;
						tmp = var_choose_order[j];
						var_choose_order[j] = var_choose_order[j + 1];
						var_choose_order[j + 1] = tmp;
					}
			/*cout << "var_implicant_count_weights_weights" << endl;
			for ( i = 0; i < core_len; i++ )
				cout << var_implicant_count_weights[i] << " ";
			cout << endl;*/

			delete[] var_implicant_count_weights;
			break;
		case 5: // PB
			for ( i = 0; i < full_mask_var_count; i++ )
				var_choose_order[i] = core_len - i - 1;
			break;
		/*case 4: // Gifford, 3 blocks, 24 vars, mask 37 (100101)
			k = 0;
			// first block 0..7
			for ( i = 0; i < Gifford_block_len; i++ )
			{
				var_choose_order[k] = i;
				k++;
			}
			// second block 16..23
			for ( i = Gifford_block_len * 2; i < Gifford_block_len * 3; i++ )
			{
				var_choose_order[k] = i;
				k++;
			}
			// third block 40..47
			for ( i = Gifford_block_len * 5; i < Gifford_block_len * 6; i++ )
			{
				var_choose_order[k] = i;
				k++;
			}
			// set -1 to unnecessary vars
			if ( k > full_mask_var_count )
				for ( i = full_mask_var_count; i < k; i++ )
					var_choose_order[i] = -1;
			else // > 24, use defaulr schema
				for ( i = 0; i < full_mask_var_count; i++ )
					var_choose_order[i] = i;
			break;*/
		case 10: // a5/1 31 = ( 11 + 11 + 9 )
			LFSR_1_arr_len = 19;
			LFSR_2_arr_len = 22;
			LFSR_3_arr_len = 23;
			LFSR_1_mid_ind = 8;
			LFSR_2_mid_ind = 29;
			LFSR_3_mid_ind = 51;
			LFSR_1_cur_ind = LFSR_1_mid_ind;
			LFSR_2_cur_ind = LFSR_2_mid_ind;
			LFSR_3_cur_ind = LFSR_3_mid_ind;
			if ( ( part_mask_var_count > 31 ) || ( part_mask_var_count == 0 ) ||
				 ( full_mask_var_count > 31 ) || ( full_mask_var_count < 4  ) )
			{
				std :: cout << "\n part_mask_var_count must be > 0 and <= 31 \n" << endl;
				return 1;
			}
			// set 4 vars of LFSR # 2 and # 3
			var_choose_order[0] = LFSR_2_cur_ind; // 27
			var_choose_order[1] = LFSR_3_cur_ind; // 51
			LFSR_2_cur_ind--;
			LFSR_3_cur_ind--;
			var_choose_order[2] = LFSR_2_cur_ind; // 26
			var_choose_order[3] = LFSR_3_cur_ind; // 50
			LFSR_2_cur_ind--;
			LFSR_3_cur_ind--;
			// check other 9 + 9 + 9 vars
			for ( i = 0; i < full_mask_var_count - 4; i++ )
			{
				switch ( i % 3 )
				{
				case 0:
					var_choose_order[i + 4] = LFSR_1_cur_ind;
					--LFSR_1_cur_ind;
					break;
				case 1:
					var_choose_order[i + 4] = LFSR_2_cur_ind;
					--LFSR_2_cur_ind;
					break;
				case 2:
					var_choose_order[i + 4] = LFSR_3_cur_ind;
					--LFSR_3_cur_ind;
					break;
				default:
					break;
				}
			}
			break;
		case 20: // first LFSR vars for sum_4_63 generator, 63 = ( 13 + 15 + 16 + 17 )
			LFSR_1_arr_len = 13;
			LFSR_2_arr_len = 15;
			LFSR_3_arr_len = 16;
			LFSR_4_arr_len = 17;
			LFSR_1_cur_ind = 0;
			LFSR_2_cur_ind = LFSR_1_arr_len;
			LFSR_3_cur_ind = LFSR_1_arr_len + LFSR_2_arr_len;
			LFSR_4_cur_ind = LFSR_1_arr_len + LFSR_2_arr_len + LFSR_3_arr_len;
			// fill var_choose_order array
			for ( i = 0; i < full_mask_var_count; i++ )
			{
				switch ( i % 4 )
				{
				case 0:
					var_choose_order[i] = LFSR_1_cur_ind;
					LFSR_1_cur_ind++;
					break;
				case 1:
					var_choose_order[i] = LFSR_2_cur_ind;
					LFSR_2_cur_ind++;
					break;
				case 2:
					var_choose_order[i] = LFSR_3_cur_ind;
					LFSR_3_cur_ind++;
					break;
				case 3:
					var_choose_order[i] = LFSR_4_cur_ind;
					LFSR_4_cur_ind++;
					break;
				default:
					break;
				}
			}
			break;
		case 21: // first LFSR vars for sum_4_66 generator, 66 = ( 13 + 15 + 17 + 19 )
			LFSR_1_arr_len = 13;
			LFSR_2_arr_len = 15;
			LFSR_3_arr_len = 17;
			LFSR_4_arr_len = 19;
			LFSR_1_cur_ind = 0;
			LFSR_2_cur_ind = LFSR_1_arr_len;
			LFSR_3_cur_ind = LFSR_1_arr_len + LFSR_2_arr_len;
			LFSR_4_cur_ind = LFSR_1_arr_len + LFSR_2_arr_len + LFSR_3_arr_len;
			// fill var_choose_order array
			for ( i = 0; i < full_mask_var_count; i++ )
			{
				switch ( i % 4 )
				{
				case 0:
					var_choose_order[i] = LFSR_1_cur_ind;
					LFSR_1_cur_ind++;
					break;
				case 1:
					var_choose_order[i] = LFSR_2_cur_ind;
					LFSR_2_cur_ind++;
					break;
				case 2:
					var_choose_order[i] = LFSR_3_cur_ind;
					LFSR_3_cur_ind++;
					break;
				case 3:
					var_choose_order[i] = LFSR_4_cur_ind;
					LFSR_4_cur_ind++;
					break;
				default:
					break;
				}
			}
			break;
		case 22: // 1st LFSR + first for other LFSR, sum_4_63 generator, 63 = ( 13 + 15 + 16 + 17 )
			LFSR_1_arr_len = 13;
			LFSR_2_arr_len = 15;
			LFSR_3_arr_len = 16;
			LFSR_4_arr_len = 17;
			LFSR_2_cur_ind = LFSR_1_arr_len;
			LFSR_3_cur_ind = LFSR_1_arr_len + LFSR_2_arr_len;
			LFSR_4_cur_ind = LFSR_1_arr_len + LFSR_2_arr_len + LFSR_3_arr_len;
			if ( full_mask_var_count < LFSR_1_arr_len ) // if var_count < 13
			{
				for ( i = 0; i < full_mask_var_count; i++ )
					var_choose_order[i] = i;
			}
			else
			{
				for ( i = 0; i < LFSR_1_arr_len; i++ )
					var_choose_order[i] = i;
				// fill var_choose_order array
				for ( i = 0; i < full_mask_var_count - LFSR_1_arr_len; i++ )
				{
					switch ( i % 3 )
					{
					case 0:
						var_choose_order[i + LFSR_1_arr_len] = LFSR_2_cur_ind;
						LFSR_2_cur_ind++;
						break;
					case 1:
						var_choose_order[i + LFSR_1_arr_len] = LFSR_3_cur_ind;
						LFSR_3_cur_ind++;
						break;
					case 2:
						var_choose_order[i + LFSR_1_arr_len] = LFSR_4_cur_ind;
						LFSR_4_cur_ind++;
						break;
					default:
						break;
					}
				}
			}
			break;
		case 30: // first LFSR vars for tresh_5_72 generator, 72 = ( 11 + 13 + 15 + 16 + 17 )
			LFSR_1_arr_len = 11;
			LFSR_2_arr_len = 13;
			LFSR_3_arr_len = 15;
			LFSR_4_arr_len = 16;
			LFSR_5_arr_len = 17;
			LFSR_1_cur_ind = 0;
			LFSR_2_cur_ind = LFSR_1_arr_len;
			LFSR_3_cur_ind = LFSR_1_arr_len + LFSR_2_arr_len;
			LFSR_4_cur_ind = LFSR_1_arr_len + LFSR_2_arr_len + LFSR_3_arr_len;
			LFSR_5_cur_ind = LFSR_1_arr_len + LFSR_2_arr_len + LFSR_3_arr_len + LFSR_4_arr_len;
			// fill var_choose_order array
			for ( i = 0; i < full_mask_var_count; i++ )
			{
				switch ( i % 5 )
				{
				case 0:
					var_choose_order[i] = LFSR_1_cur_ind;
					LFSR_1_cur_ind++;
					break;
				case 1:
					var_choose_order[i] = LFSR_2_cur_ind;
					LFSR_2_cur_ind++;
					break;
				case 2:
					var_choose_order[i] = LFSR_3_cur_ind;
					LFSR_3_cur_ind++;
					break;
				case 3:
					var_choose_order[i] = LFSR_4_cur_ind;
					LFSR_4_cur_ind++;
					break;
				case 4:
					var_choose_order[i] = LFSR_5_cur_ind;
					LFSR_5_cur_ind++;
					break;
				default:
					break;
				}
			}
			break;
		case 31: // first LFSR vars for tresh_5_80 generator, 80 = ( 13 + 15 + 16 + 17 + 19 )
			LFSR_1_arr_len = 13;
			LFSR_2_arr_len = 15;
			LFSR_3_arr_len = 16;
			LFSR_4_arr_len = 17;
			LFSR_5_arr_len = 19;
			LFSR_1_cur_ind = 0;
			LFSR_2_cur_ind = LFSR_1_arr_len;
			LFSR_3_cur_ind = LFSR_1_arr_len + LFSR_2_arr_len;
			LFSR_4_cur_ind = LFSR_1_arr_len + LFSR_2_arr_len + LFSR_3_arr_len;
			LFSR_5_cur_ind = LFSR_1_arr_len + LFSR_2_arr_len + LFSR_3_arr_len + LFSR_4_arr_len;
			// fill var_choose_order array
			for ( i = 0; i < full_mask_var_count; i++ )
			{
				switch ( i % 5 )
				{
				case 0:
					var_choose_order[i] = LFSR_1_cur_ind;
					LFSR_1_cur_ind++;
					break;
				case 1:
					var_choose_order[i] = LFSR_2_cur_ind;
					LFSR_2_cur_ind++;
					break;
				case 2:
					var_choose_order[i] = LFSR_3_cur_ind;
					LFSR_3_cur_ind++;
					break;
				case 3:
					var_choose_order[i] = LFSR_4_cur_ind;
					LFSR_4_cur_ind++;
					break;
				case 4:
					var_choose_order[i] = LFSR_5_cur_ind;
					LFSR_5_cur_ind++;
					break;
				default:
					break;
				}
			}
			break;
		case 32: // 1st LFSR + first vars for other LFSR, tresh_5_72 generator, 72 = (11+13+15+16+17)
			LFSR_1_arr_len = 11;
			LFSR_2_arr_len = 13;
			LFSR_3_arr_len = 15;
			LFSR_4_arr_len = 16;
			LFSR_5_arr_len = 17;
			LFSR_2_cur_ind = LFSR_1_arr_len;
			LFSR_3_cur_ind = LFSR_1_arr_len + LFSR_2_arr_len;
			LFSR_4_cur_ind = LFSR_1_arr_len + LFSR_2_arr_len + LFSR_3_arr_len;
			LFSR_5_cur_ind = LFSR_1_arr_len + LFSR_2_arr_len + LFSR_3_arr_len + LFSR_4_arr_len;
			// fill var_choose_order array
			if ( full_mask_var_count < LFSR_1_arr_len ) // if var_count < 13
			{
				for ( i = 0; i < full_mask_var_count; i++ )
					var_choose_order[i] = i;
			}
			else
			{
				for ( i = 0; i < LFSR_1_arr_len; i++ )
					var_choose_order[i] = i;
				// fill var_choose_order array
				for ( i = 0; i < full_mask_var_count - LFSR_1_arr_len; i++ )
				{
					switch ( i % 4 )
					{
					case 0:
						var_choose_order[i + LFSR_1_arr_len] = LFSR_2_cur_ind;
						LFSR_2_cur_ind++;
						break;
					case 1:
						var_choose_order[i + LFSR_1_arr_len] = LFSR_3_cur_ind;
						LFSR_3_cur_ind++;
						break;
					case 2:
						var_choose_order[i + LFSR_1_arr_len] = LFSR_4_cur_ind;
						LFSR_4_cur_ind++;
						break;
					case 3:
						var_choose_order[i + LFSR_1_arr_len] = LFSR_5_cur_ind;
						LFSR_5_cur_ind++;
						break;
					default:
						break;
					}
				}
			}
			break;
		default:
			printf( "\n Unknown value of schema_type, changed to 0" );
			for ( i = 0; i < full_mask_var_count; i++ )
				var_choose_order[i] = i;
			break;
	}

	cout << "var_choose_order" << endl;
	for ( i = 0; i < core_len; i++ )
		cout << var_choose_order[i] << " ";

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
	string line_str = "", word_str = "";
	bool IsUncorrectLine = false;
	int lit_val = 0, sign = 0, lit_index = 0, val = 0;
	
	// check file with main CNF
	ifstream main_cnf( input_cnf_name, ios::in );
    if ( !main_cnf )
	{
		std :: cerr << endl << "Error in opening of file with input CNF " 
			        << input_cnf_name << endl;
		return false;
	}

	// step 1 - get var_count and clause_count;
    while ( getline( main_cnf, line_str ) )
	{
		if ( ( line_str[0] == 'p' ) || ( line_str[0] == 'c' ) )
			continue;

		else // try to read line with clause 
		{
			current_lit_count = 0;
			line_str = " " + line_str; // add space to line for correct work of parser
			for ( i = 0; i < line_str.length( ) - 1; i++ )
			{
				IsUncorrectLine = false;
				if ( ( line_str[i] == ' ' ) && ( line_str[i + 1] != ' ' ) && 
					 ( line_str[i + 1] != '0' ) )
				{
					word_str = ""; // init string for cuttenr word from string
					k = i + 1; // k = index of first symbol in current word
					do
					{
						word_str += line_str[k];
						k++;
						if ( k == line_str.length( ) ) // skip empty or uncorrect line
						{
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
	lits_clause_lengths = new int[lit_count];
	for ( i = 0; i < ( unsigned int )lit_count; i++ )
		lits_clause_lengths[i] = 0;

	main_cnf.close( ); // reopen file
	main_cnf.clear( );
	main_cnf.open( input_cnf_name );
	current_clause_count = 0;

	// step 2 - get arrays of lengths
	while ( getline( main_cnf, line_str ) )
	{
		if ( ( line_str[0] == 'p' ) || ( line_str[0] == 'c' ) )
			continue;

		else // try to read line with clause 
		{
			current_lit_count = 0;
			line_str = " " + line_str; // add space to line for correct work of parser
			for ( i = 0; i < line_str.length( ) - 1; i++ )
			{
				IsUncorrectLine = false;
				if ( ( line_str[i] == ' ' ) && ( line_str[i + 1] != ' ' ) && 
					 ( line_str[i + 1] != '0' ) )
				{
					word_str = ""; // init string for cuttenr word from string
					k = i + 1; // k = index of first symbol in current word
					do
					{
						word_str += line_str[k];
						k++;
						if ( k == line_str.length( ) ) // skip empty or uncorrect line
						{
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
					if ( lit_val < 0 )
					{
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
	int lit_val = 0,
		i = 0,
		input_var_num = -1,
		first_obj_var = -1;
	string line_str = "", 
		   word_str = "";
	bool IncorrectLine = false;

	if ( !ReadVarCount( ) ) {
		cerr << "Error in ReadVarCount" << endl;
		return false;
	}
	//cout << "*** Success of ReadVarCount" << endl;

	clause_array      = new int*[clause_count];
	lits_clause_array = new int*[lit_count]; // array of indexes of clauses for literals

	int *lits_clause_current = new int[lit_count];

	for ( i = 0; i < clause_count; i++ )
		clause_array[i]      = new int[clause_lengths[i]];
	for ( i = 0; i < lit_count; i++ )
	{
		lits_clause_array[i] = new int[lits_clause_lengths[i]];
		lits_clause_current[i] = 0;
	}

	//std :: cout << "\n***3" << endl;
	// check file with main CNF
	ifstream main_cnf( input_cnf_name, ios::in );
    if ( !main_cnf )
	{
		cerr << "Error in opening of file with input CNF with name" 
			 << input_cnf_name << endl;
		return false;
	}

	current_clause_count = 0;
	// if core_len wasn't set by key parameter then by default core_len = count of all variables
	if ( !core_len ) {
		if ( var_count < MAX_CORE_LEN ) core_len = var_count;
		else core_len = MAX_CORE_LEN;
	}

	stringstream sstream;
	string str1 = "", str2 = "", str3 = "", str4 = "";
	bool Is_InpVar = false, Is_ConstrLen = false, Is_ObjLen = false, Is_ObjVars = false;
	while ( getline( main_cnf, line_str ) )
	{
		if ( line_str[0] == 'p' )
			continue;

		if ( line_str[0] == 'c' ) // in comment string can exist count of input variables
		{
			//parse string for ex. "c 1452 input variables" or "c input variables 1452"
			sstream << line_str;
			sstream >> str1 >> str2 >> str3 >> str4; // get and parse words in string
			sstream.str( "" );
			sstream.clear( );

			if ( !Is_InpVar )
			{
				if ( ( str2 == "input" ) && ( str3 == "variables" ) ) 
					istringstream( str4 ) >> input_var_num;
				else if ( ( str3 == "input" ) && ( str4 == "variables" ) ) 
					istringstream( str2 ) >> input_var_num;
				if ( input_var_num > 0 )
				{
					core_len = input_var_num;
					if ( core_len > MAX_CORE_LEN ) {
						core_len = MAX_CORE_LEN;
						cout << "Warning. core_len > MAX_CORE_LEN " << MAX_CORE_LEN << " Changed to MAX_CORE_LEN";
					}
					Is_InpVar = true;
					continue;
				}
			}
			if ( !Is_ConstrLen )
			{
				if ( ( str2 == "constraint" ) && ( str3 == "clauses" ) ) 
				{
					istringstream( str4 ) >> constr_clauses_count; 
					Is_ConstrLen = true;
					continue;
				}
			}
			if ( !Is_ObjLen )
			{
				if ( ( ( str2 == "object" ) || ( str2 == "obj" ) ) && ( str3 == "clauses" )  )
				{
					istringstream( str4 ) >> obj_clauses_count;
					Is_ObjLen = true;
					continue;
				}
			}
			if ( !Is_ObjVars )
			{
				if ( ( str2 == "obj" ) && ( str3 == "vars" ) ) 
					istringstream( str4 ) >> first_obj_var;
				else if ( ( str2 == "object" ) && ( str3 == "variables" ) ) 
					istringstream( str4 ) >> first_obj_var;
				if ( first_obj_var > 0 )
				{
					if ( !Is_ConstrLen ) constr_clauses_count = clause_count;

					// read obj vars indexes from string "c obj vars ...")
					obj_vars_count = 0;
					line_str += ' ';
					while ( line_str.length( ) )
					{
						while ( ( line_str.length( ) > 0 ) && ( line_str[0] == ' ' ) ) 
							line_str.erase( 0, 1 ); // skip spaces
						if ( !line_str.length( ) )
							break;
						int space_pos = line_str.find( " " ); // find space after word
						word_str = line_str.substr( 0, space_pos ); // get word
						line_str.erase( 0, space_pos + 1 ); // delete word and space after it
						int word_value = atoi( word_str.c_str( ) );
						if ( word_value > 0 ) // if it is number > 0
						{
							if ( obj_vars_count < 32 )
							{
								obj_vars[obj_vars_count] = word_value; // add to array
								obj_vars_count++;
							}
							else
								cout << "\n ***Error. String c obj vars ... contains too many values ";
						}
					}
					Is_ObjVars = true;
					continue;
				} // if ( sscanf( line_str.c_str( ), "c obj vars %d", &first_obj_var ) )
			}  // if ( !Is_ObjVars )
		}
		else // if ( ( line_str[0] != 'p' ) && ( line_str[0] != 'c' ) )
		{
			// try to read line with clause
			current_lit_count = 0; // cuttenr count of lits in current clause
			line_str = " " + line_str;
			for ( i = 0; i < ( int )line_str.length( ) - 1; i++ )
			{
				IncorrectLine = false;
				if ( ( line_str[i] == ' ' ) && ( line_str[i + 1] != ' ' ) && 
					 ( line_str[i + 1] != '0' ) )
				{
					word_str = "";
					k = i + 1;
					do
					{
						word_str += line_str[k];
						k++;
						if ( k == line_str.length( ) ) // skip empty or uncorrect line
						{
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
	
	// create array for answer
	b_SAT_set_array = new int[var_count];
	for ( i = 0; i < var_count; i++ )
		b_SAT_set_array[i] = -1;

	// if PB data is correct then turn on PB mode
	if ( ( constr_clauses_count > 0 ) && ( obj_vars_count > 0 ) ) 
	{
		IsPB = true;
		solver_type = 3;
		/*if ( ( best_lower_bound > -1 ) && ( upper_bound > 0 ) ) // if bounds then equality mode
		{
			PB_mode = 2;
			if ( best_lower_bound > upper_bound )  // check correctness
				best_lower_bound = upper_bound;
		}
		else PB_mode = 1;*/
	}

	delete[] lits_clause_current;
	main_cnf.close( );

	return true;
}

//---------------------------------------------------------
bool MPI_Base :: CheckSATset( int &int_answer )
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
	while ( ( cnf_sat_val ) && ( checked_clauses != clause_count ) )
	{
		// add one to count of checked clauses in main CNF
        checked_clauses++;
		clause_sat_val = 0;
		j = 0;
		// check literals while all of checked is UNSAT
		while ( ( !clause_sat_val ) && ( j < clause_lengths[i] ) )
		{
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
	// take bool answer
	int_answer = cnf_sat_val;
	// if OK
	return true;
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
		i, k,
		answer_var_count;
	bool bIsSATSetExist = false;
	string answer_file_name,
		   output_file_name = "output",
		   str_answer,
		   line_buffer;
	ofstream answer_file,
			 output_file;

	cout << "\n Start of AnalyzeSATset";
	lit_SAT_set_array = new int[var_count];
	for ( i = 0; i < ( int )var_count; i++ )
	{
		sign = 0;
		if ( !b_SAT_set_array[i] ) 
			sign = 1;
		// ( b_SAT_set_array[0] = 1 ) -> ( lit_SAT_set_array[0] = 2 )
		val = ( ( i + 1 ) << 1 ) + sign; // ( b_SAT_set_array[0] = 0 ) -> ( lit_SAT_set_array[0] = 3 )
		lit_SAT_set_array[i] = val;
	}
	//std :: cout << "\n***Start CheckSATset" << endl;
	// if SAT set exist then check it
	if ( !CheckSATset( int_answer ) )
	{
		cerr << "Error in checking of SAT set" << endl;
		// deallocate arrays
		/*
		delete[] clause_array;
		delete[] clause_lengths;
		*/
		//
		return false;
	}
	// 
	if ( int_answer ) {
		str_answer = "\n SAT set is correct";
		IsSAT = true;
	}
	else {
		str_answer = "\n SAT set is incorrect ";
		IsSAT = false;
	}
	cout << endl << "str_answer is " << str_answer;

	// open file for writing of answer in zchaff style
	//str_answer = "answer_";
	answer_file_name = "answer_";
	answer_file_name += input_cnf_name;
	
	k = 1;
	for( i = 1; i < ( int )answer_file_name.length( ); i++ )
	{
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
	answer_file.open( answer_file_name.c_str( ), ios :: out );
	if ( !( answer_file.is_open( ) ) )
	{
		std :: cerr << "Error in opening of file with answer in zchaff style " 
			        << answer_file_name << endl;
		return false;
	}
	// open file for writing of output info
	output_file.open( output_file_name.c_str( ), ios :: out );
	output_file.write( str_answer.c_str( ), str_answer.length( ) );
	// wtite answer to file
	stringstream sstream;
	sstream << endl;
	for ( i = 0; i < ( int )var_count; i++ )
		sstream << "x" << i + 1 << "=" << b_SAT_set_array[i] << " ";

	output_file << sstream.rdbuf( );

	// put only core vars if such was specified
	if ( !core_len ) answer_var_count = var_count;
	else answer_var_count = core_len;

	sstream.str( "" ); // clear stringstream
	sstream.clear();
	for ( i = 0; i < answer_var_count; i++ )
	{
		if ( b_SAT_set_array[i] ) sstream << i + 1;
		else sstream << -( i + 1 );
		sstream << " ";
	}
	answer_file << sstream.rdbuf( );
	// close file
	output_file.close( );
	answer_file.close( );

	// deallocate memory
	delete[] lit_SAT_set_array;
	// if OK
	return true;
}

//---------------------------------------------------------
bool MPI_Base :: cpuTimeInHours( double full_seconds, int &real_hours, int &real_minutes, int &real_seconds ) 
{
// Time of work in hours, minutes, seconds
	int full_minutes = -1;
	//
	full_minutes = ( int )full_seconds / 60;
	real_seconds = ( int )full_seconds % 60;
	real_hours   = full_minutes / 60;
	real_minutes = full_minutes % 60;
	//
	return true;
}

//---------------------------------------------------------
bool MPI_Base :: WriteTimeToFile( double whole_time_sec )
{
// Write time of solving to output file
	string output_file_name = "output",
		   str_output,
		   line_buffer,
		   TimeIs_char = "\n Time: ",
		   new_line_char = "\n",
		   hour_char = " h ",
		   minute_char = " m ",
		   second_char = " s ";
	int str_output_cur_ind = 0,
		real_hours = -1,
		real_minutes = -1,
		real_seconds = -1;
	ofstream output_file;
	//
	// std :: cout << "\n time in sec " << whole_time_sec << endl;
	// open file for writing of output info
	output_file.open( output_file_name.c_str( ), ios :: app );
	if ( !output_file.is_open( ) ) {
		std :: cout << "Error in opening of output file " << output_file_name <<  endl;
		return false;
	}
	//
	if ( !cpuTimeInHours( whole_time_sec, real_hours, real_minutes, real_seconds ) ) 
		{ cout << "Error in cpuTimeInHours" << endl; return false; }

	if ( IsSAT ) str_output = " SAT \n";
	else str_output = " UNSAT \n";

	stringstream sstream;
	sstream << str_output << TimeIs_char << real_hours << hour_char 
		    << real_minutes << minute_char << real_seconds << second_char;
	output_file << sstream.rdbuf( );
	
	output_file.close( );

	return true;
}

//---------------------------------------------------------
bool MPI_Base :: SortValuesDecrease( unsigned int range, 
									 unsigned int *&sorted_index_array )
{
// Sorting of array with values
// by decreasing of weight (at first - 1-vector, at last 0-vector) 
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

	for ( i = 0; i < bit_value_count; i++ ) {
		//index_array[i] = new unsigned int[MAX_INDEX_ARRAY_LENGTH];
		index_array_lengths[i] = 0;
	}
	
	// compute lengths of arrays, they are same for any kind of values
	for ( i = 0; i < table_val_count; i++ )
	{
		// get count of "1" in 32-bit value
		// shift 3 times
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
		// null array of lengths
		index_array_lengths[i] = 0;
	}
	
	// printf( "\n table_val_count is %d", table_val_count );
	// fill array of indexes
	for ( i = 0; i < table_val_count; i++ )
	{
		// by shifting 3 times get count of "1" in 32-bit value
		bit_weight = bits[( unsigned char )( i )] + bits[( unsigned char )( i >> 8 )] +
		             bits[( unsigned char )( i >> 16 )] + bits[( unsigned char )( i >> 24 )];
		// again fill array of lengths
		index_array_lengths[bit_weight] += 1;
		index_array[bit_weight][index_array_lengths[bit_weight] - 1] = i;
		// printf( "\n done # %d", i );
	}

	k = 0;
	// fill array with sorted values
	for ( i = bit_value_count - 1; i > -1; --i )
	{
		val2 = index_array_lengths[i];
		for ( j = 0; j < val2; j++ ) {
			sorted_index_array[k] = index_array[i][j];
			k++;
		}
	}
	
	// deallocate arrays
	free( sort_index_arr );

	for ( i = 0; i < bit_value_count; i++ )
	{
		if ( index_array_lengths[i] )
			free( index_array[i] );
	}
	
	free( index_array_lengths );
	free( index_array );
	
	return true;
}