#ifdef _MPI
#include <mpi.h>
#endif

#ifdef _WIN32
#define inline __inline // compatible with MS VS 6.0
#endif

#include "./dminisat/dminisat_vec.h"

//=================================================================================================
// Simple types:

typedef int                lit;
typedef char               lbool_dminisat;

static const lbool_dminisat l_Undef_dminisat   =  0;
static const lbool_dminisat l_True_dminisat    =  1;
static const lbool_dminisat l_False_dminisat   = -1;

// does not work for c++
typedef int boolean_dminisat;
static const boolean_dminisat ltrue = 1;
static const boolean_dminisat lfalse = 0;

#ifdef _WIN32
typedef signed __int64     uint64;   // compatible with MS VS 6.0
#else
typedef unsigned long long uint64;
#endif

static const int   var_Undef = -1;
static const lit   lit_Undef_dminisat = -2;

static inline lit  toLit_dminisat   (int v) { return v + v; }
static inline lit  lit_neg (lit l) { return l ^ 1; }
static inline int  lit_var (lit l) { return l >> 1; }
static inline int  lit_sign(lit l) { return (l & 1); }

//=================================================================================================
// Public interface:

struct solver_t;
typedef struct solver_t solver;

extern solver* solver_new(void);
extern void    solver_delete(solver* s);

extern boolean_dminisat    solver_addclause(solver* s, lit* begin, lit* end);
extern boolean_dminisat    solver_simplify(solver* s);
extern boolean_dminisat    solver_solve( solver* s, lit* begin, lit* end );

extern int     solver_nvars( solver* s );
extern int     solver_nclauses( solver* s );
extern int     solver_nconflicts( solver* s );

extern void    solver_setnvars( solver* s, int n );

extern void solver_cleardb( solver* s );

struct stats_t
{
    uint64   starts, decisions, propagations, inspects, conflicts;
    uint64   clauses, clauses_literals, learnts, learnts_literals, max_literals, tot_literals;
};
typedef struct stats_t lstats;

//=================================================================================================
// Solver representation:

struct clause_t;
typedef struct clause_t clause;

struct solver_t
{
    // added paremeters
	int IsPredict;
	int core_len;
	double corevars_activ_type;
	int IsHardProblem;

	// native parameters
	int      size;          // nof variables
    int      cap;           // size of varmaps
    int      qhead;         // Head index of queue.
    int      qtail;         // Tail index of queue.

    // clauses
    vecp     clauses;       // List of problem constraints. (contains: clause*)
    vecp     learnts;       // List of learnt clauses. (contains: clause*)

    // activities
    double   var_inc;       // Amount to bump next variable with.
    double   var_decay;     // INVERSE decay factor for variable activity: stores 1/decay.
    float    cla_inc;       // Amount to bump next clause with.
    float    cla_decay;     // INVERSE decay factor for clause activity: stores 1/decay.

    vecp*    wlists;        //
    double*  activity;      // A heuristic measurement of the activity of a variable.
    lbool_dminisat*   assigns;       // Current values of variables.
    int*     orderpos;      // Index in variable order.
    clause** reasons;       //
    int*     levels;        //
    lit*     trail;

    clause*  binary;        // A temporary binary clause
    lbool_dminisat*   tags;          //
    veci     tagged;        // (contains: var)
    veci     stack;         // (contains: var)

    veci     order;         // Variable order. (heap) (contains: var)
    veci     trail_lim;     // Separator indices for different decision levels in 'trail'. (contains: int)
    veci     model;         // If problem is solved, this vector contains the model (contains: lbool_dminisat).

    int      root_level;    // Level of first proper decision.
    int      simpdb_assigns;// Number of top-level assignments at last 'simplifyDB()'.
    int      simpdb_props;  // Number of propagations before next 'simplifyDB()'.
    double   random_seed;
    double   progress_estimate;
    int      verbosity;     // Verbosity level. 0=silent, 1=some progress report, 2=everything

    lstats    stats;
};

