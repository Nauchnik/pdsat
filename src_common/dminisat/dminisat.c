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
    char*   data = (char*)malloc( 65536 );
    int     cap  = 65536;
    int     size = 0;

    while ( !feof( in ) )
	{
        if ( size == cap )
		{
            cap *= 2;
            data = (char*)realloc( data, cap ); 
		}
        size += fread( &data[size], 1, 65536, in );
    }
    data = (char*)realloc( data, size + 1 );
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
