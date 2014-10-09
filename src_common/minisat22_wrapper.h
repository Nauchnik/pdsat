#ifndef minisat22_wrapper_h
#define minisat22_wrapper_h

#include <errno.h>

#include <signal.h>
#include <iostream>
#include <vector>

#include "utils/System.h"
#include "utils/ParseUtils.h"
#include "utils/Options.h"
#include "core/Dimacs.h"
#include "core/Solver.h"
//#include "simp/SimpSolver.h"

using namespace Minisat;

class minisat22_wrapper
{
public:
	void readClause(StreamBuffer& in, vec<Lit>& lits);
	void parse_DIMACS_to_problem(std::istream& input, Problem& cnf);
	void parse_DIMACS_from_inc( std::vector<int> &cnf_vec, Problem& cnf);
	void printProblem(const Problem& p, std::ostream& out);
};

#endif