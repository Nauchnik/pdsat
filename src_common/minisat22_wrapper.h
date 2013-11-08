#ifndef minisat22_wrapper_h
#define minisat22_wrapper_h

#include <errno.h>

#include <signal.h>
#include <iostream>
#include <fstream>
#include <vector>

#include "minisat/utils/System.h"
#include "minisat/utils/ParseUtils.h"
#include "minisat/utils/Options.h"
#include "minisat/core/Dimacs.h"
#include "minisat/core/Solver.h"
#include "minisat/simp/SimpSolver.h"

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