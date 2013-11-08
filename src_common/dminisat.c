// dminisat

#include <stdio.h>
#include <time.h>

#include "dminisat.h"

#pragma warning( disable : 4996 )

int proc_status;
time_t start_time;

//=================================================================================================
// Helpers:


// Reads an input stream to end-of-file and returns the result as a 'char*' terminated by '\0'
// (dynamic allocation in case 'in' is standard input).
//

static char* readFile( FILE *  in )
{
    char*   data = malloc( 65536 );
    int     cap  = 65536;
    int     size = 0;

    while ( !feof( in ) )
	{
        if ( size == cap )
		{
            cap *= 2;
            data = realloc( data, cap ); 
		}
        size += fread( &data[size], 1, 65536, in );
    }
    data = realloc( data, size + 1 );
    data[size] = '\0';
	
    return data;
}

//static inline double cpuTime(void) {
//    struct rusage ru;
//    getrusage(RUSAGE_SELF, &ru);
//    return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000; }

//=================================================================================================
// DIMACS Parser:

static inline void skipWhitespace(char** in) 
{
    while ((**in >= 9 && **in <= 13) || **in == 32)
        (*in)++; 
}

static inline void skipLine(char** in) 
{
    for (;;)
	{
        if (**in == 0) return;
        if (**in == '\n') { (*in)++; return; }
        (*in)++; 
	} 
}

static inline int parseInt(char** in) 
{
    int     val = 0;
    int    _neg = 0;
    skipWhitespace(in);
    if      (**in == '-') _neg = 1, (*in)++;
    else if (**in == '+') (*in)++;
    if (**in < '0' || **in > '9') fprintf(stderr, "PARSE ERROR! Unexpected char: %c\n", **in), exit(1);
    while (**in >= '0' && **in <= '9')
        val = val*10 + (**in - '0'),
        (*in)++;
    return _neg ? -val : val; 
}

static void readClause(char** in, solver* s, veci* lits) 
{
    int parsed_lit, var;
    veci_resize(lits,0);
    for (;;)
	{
        parsed_lit = parseInt(in);
        if (parsed_lit == 0) break;
        var = abs(parsed_lit)-1;
        veci_push(lits, (parsed_lit > 0 ? toLit_dminisat(var) : lit_neg(toLit_dminisat(var))));
    }
}

static lbool_dminisat parse_DIMACS_main(char* in, solver* s) 
{
    veci lits;
    veci_new(&lits);

    for ( ; ; )
	{
        skipWhitespace( &in );
        if (*in == 0)
            break;
/*
        else if (*in == 'p') {
			int nvars = 0;

			sscanf(in, "p cnf %d", &nvars);
			solver_setnvars(s, nvars + 1);
            skipLine(&in);
		}
*/
        else if (*in == 'c' || *in == 'p')
            skipLine(&in);
        else
		{
            lit* begin;
            readClause(&in, s, &lits);
            begin = veci_begin(&lits);
            if (!solver_addclause(s, begin, begin+veci_size(&lits)))
			{
                veci_delete(&lits);
                return l_False_dminisat;
            }
        }
    }
    veci_delete(&lits);
    return solver_simplify(s);
}


// Inserts problem into solver. Returns FALSE upon immediate conflict.
//
lbool_dminisat parse_DIMACS( FILE * in, solver* s ) 
{
    char* text = readFile( in );
    lbool_dminisat ret  = parse_DIMACS_main( text, s );
    free( text );
    return ret; 
}

//=================================================================================================

static void printStats(lstats* stats, int cpu_time)
{
    double Time    = (float)(cpu_time)/(float)(CLOCKS_PER_SEC);
    printf("restarts          : %12d\n", stats->starts);
    printf("conflicts         : %12.0f           (%9.0f / sec      )\n",  (double)stats->conflicts   , (double)stats->conflicts   /Time);
    printf("decisions         : %12.0f           (%9.0f / sec      )\n",  (double)stats->decisions   , (double)stats->decisions   /Time);
    printf("propagations      : %12.0f           (%9.0f / sec      )\n",  (double)stats->propagations, (double)stats->propagations/Time);
    printf("inspects          : %12.0f           (%9.0f / sec      )\n",  (double)stats->inspects    , (double)stats->inspects    /Time);
    printf("conflict literals : %12.0f           (%9.2f %% deleted  )\n", (double)stats->tot_literals, (double)(stats->max_literals - stats->tot_literals) * 100.0 / (double)stats->max_literals);
    printf("CPU time          : %12.2f sec\n", Time);
}

//solver* slv;
//static void SIGINT_handler(int signum) {
//    printf("\n"); printf("*** INTERRUPTED ***\n");
//    printStats(&slv->stats, cpuTime());
//    printf("\n"); printf("*** INTERRUPTED ***\n");
//    exit(0); }
/*
//=================================================================================================
static unsigned count_ones( long long num )
{
	unsigned cnt = 0;
	
	int i;
	for ( i = 0; i < NUM_KEY_BITS; i++ ) 
	{
		if ( num & 1 )
			cnt++;

		num >>= 1;
	}

	return cnt;
}
*/
//=================================================================================================
static unsigned new_count_ones( unsigned num )
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

//---------------------------------------------------------
static unsigned int SortValues( unsigned int range, 
							    unsigned int *sorted_index_array )
{
// Sorting of array with values
// from most probably (with equal number of 0 and 1) 
// to in fact inpossible values (0- and 1-vector)
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
	unsigned int table_val_count = 0,
				 i = 0,
				 j = 0,
				 bit_weight = 0,
				 bit_value_count = 0,
				 mid_val = 0,
				 k,
				 val1, 
				 val2, 
				 val3;
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
	
	// compute lengths of arrrays, they are same for any kind of values
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
	mid_val = bit_value_count / 2;

	// fill array with all indexes sorted by decreasing of probability
	if ( bit_value_count % 2 )
	{
		// example: 0 1 2 -> 1 0 2
		sort_index_arr[k] = mid_val;
		k++;
		for ( i = 0; i < mid_val; i++ )
		{
			sort_index_arr[k]     = mid_val - i - 1;
			sort_index_arr[k + 1] = mid_val + i + 1;
			k = k + 2;
		}
	}	
	else
	{
		// example: 0 1 2 3 -> 1 2 0 3
		for ( i = 0; i < mid_val; i++ )
		{
			sort_index_arr[k]     = mid_val - i - 1;
			sort_index_arr[k + 1] = mid_val + i;
			k = k + 2;
		}
	}

	k = 0;
	// fill array with sorted values
	for ( i = 0; i < bit_value_count; i++ )
	{
		val1 = sort_index_arr[i];
		val2 = index_array_lengths[val1];
		for ( j = 0; j < val2; j++ )
		{
			val3 = index_array[val1][j];
			sorted_index_array[k] = val3;
			val3 = sorted_index_array[k];
			k++;
		}
	}
	
	// deallocate arrays
	free( sort_index_arr );

	for ( i = 0; i < bit_value_count; i++ )
		if ( index_array_lengths[i] )
			free( index_array[i] );
	
	free( index_array_lengths );
	free( index_array );
	
	return 1;
}

//=================================================================================================
void shl64( unsigned long long int *val_for_left_shift, unsigned int bit_count )
{
	unsigned int val1, val2; 
	if ( ( bit_count > 31 ) && ( bit_count <= 62 ) )
	{
		val1 = 31;
		val2 = bit_count - val1;
		*val_for_left_shift =  ( unsigned long long int )( 1 << val1 );
		*val_for_left_shift *= ( unsigned long long int )( 1 << val2 );
	}
	else if ( bit_count <= 31 )
		*val_for_left_shift =  ( unsigned long long int )( 1 << bit_count );
}

//---------------------------------------------------------
int dminisat_solve( char* CNFfile, unsigned int full_mask[FULL_MASK_LEN], 
			 unsigned int part_mask[FULL_MASK_LEN], unsigned int values[FULL_MASK_LEN], 
			 int solver_type, int core_len, double start_activity, int *process_sat_count, 
			 int **b_SAT_set_array, int sort_type, double *cnf_time_from_node,
			 int IsPredict, int IsHardProblem )
{
// Solve SAT-problem with solver dminisat
    solver *s = solver_new( );
    lbool_dminisat st;
    FILE *in;
    lit assigns[MAX_ASSIGNS_COUNT];
    // full_mask; // полная маска, единичными битами помечены литералы, которые передаются в solver_solve как известные
    // part_mask; // частичная маска, единичные биты соответствуют литералам, которые переданы солверу как известные
    // values; // значения литералов, указанных в part_mask
	// full_mask_2; // полная маска # 2
    // part_mask_2; // частичная маска # 2
    // values_2; // значения литералов, указанных в part_mask # 2
	unsigned long long int range_val_count,
				  lint,
				  range_mask;
	unsigned int  range = 0,
				  nassigns,
				  i, j,
				  clk,
				  mask;
	unsigned int range_mask_ind;
	unsigned int cur_var_ind;

	clk = clock( );
    in = fopen( CNFfile, "rb" );
	
	if ( in == NULL ) {
		printf( "\nERROR! Could not open file: %s", CNFfile );
		return 0;
	}

    st = parse_DIMACS( in, s );
    fclose( in );

    if ( st == l_False_dminisat ) {
        solver_delete( s );
        printf( "\nTrivial problem\nUNSATISFIABLE" );
        return 0;
    }
	
    s->verbosity = 0;
	
	for ( i = 1; i < FULL_MASK_LEN; i++ )
		range += new_count_ones( full_mask[i] ^ part_mask[i] );

	shl64( &range_val_count, range );

	for ( lint = 0; lint < range_val_count; lint++ ) {
		nassigns = 0;
		range_mask = 1;
		range_mask_ind = 1;

		for ( i = 1; i < FULL_MASK_LEN; i++ ) {
			for ( j = 0; j < UINT_LEN; j++ ) {
				mask = ( 1 << j );
				if ( part_mask[i] & mask ) {
					cur_var_ind = ( i-1 ) * UINT_LEN + j;
					if ( values[i] & mask )
						assigns[nassigns] = toLit_dminisat( cur_var_ind );
					else
						assigns[nassigns] = lit_neg( toLit_dminisat( cur_var_ind ) );

					nassigns++;
				}
				else if ( full_mask[i] & mask ) {
					cur_var_ind = ( i-1 ) * UINT_LEN + j;
					if ( lint & range_mask )
						assigns[nassigns] = toLit_dminisat( cur_var_ind  );
					else
						assigns[nassigns] = lit_neg( toLit_dminisat( cur_var_ind  ) );

 					shl64( &range_mask, range_mask_ind );
					
					range_mask_ind++;
					nassigns++;
				}
			}
		}
	
		/*printf( "\n nassigns %d", nassigns );
		for ( i = 0; i < nassigns; i++ ) 
			printf( "\n assigns %d", assigns[i] );*/

		//printf( "\n Start of solver_solve" );
		s->IsPredict           = IsPredict;
		s->start_activity      = start_activity;
		s->core_len			   = core_len;
		s->IsHardProblem       = IsHardProblem;
		s->verbosity		   = 0; // no info

#ifdef _MPI
	if ( IsPredict ) // if predicting then get start time-label 
		*cnf_time_from_node = MPI_Wtime( ); // get first time of solving
#endif
		st = solver_solve( s, assigns, assigns + nassigns );

#ifdef _MPI
	if ( IsPredict ) // if predicting then get real time
		*cnf_time_from_node = MPI_Wtime( ) - *cnf_time_from_node;
#endif
		//printf( "\n End of solver_solve" );
		//printf( "\n 4" );
		// clean learnts clauses to start again
		solver_cleardb( s );
		
		if ( st == l_True_dminisat )
	    	break; // exit form bottom for( )
			
	} // for ( lint = 0; lint < range_val_count; lint++ ) 

	// print set assignment
	if ( st == l_True_dminisat ) {
		int k;
		printStats( &s->stats, clock( ) - clk );
		printf( "Satisfying solution: \n" );
		*(process_sat_count)++;
		for ( k = 0; k < s->model.size; k++ ) {
			printf( "x%d=%d ", k, s->model.ptr[k] == l_True_dminisat );
			if ( s->model.ptr[k] == l_True_dminisat )
				b_SAT_set_array[0][k] = 1;
			else
				b_SAT_set_array[0][k] = 0;
		}
	}
	
    solver_delete( s );

    return 1;
}
