#ifndef MINISAT_LOGGER_H
#define MINISAT_LOGGER_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include "minisat/mtl/Vec.h"
#include "minisat/core/SolverTypes.h"

void WriteLearntClauses(const Minisat::vec<Minisat::CRef>& learnts, const Minisat::ClauseAllocator& ca)
{
	static int file_id = 0;
	std::stringstream ss;
	ss << "LearntsClauses_" << ++file_id << ".cnf";
	std::ofstream out(ss.str().c_str(), std::ios::out);
	if(out.is_open())
	{
		for(int i = 0; i < learnts.size(); i++)
		{
			const Minisat::Clause& c = ca[learnts[i]];
			for(int j = 0; j < c.size(); j++)
			{
				const Minisat::Lit& lit = c[j];
				out << (Minisat::sign(lit) ? "-" : "") << Minisat::var(lit) << " ";
			}
			out << "0\n";
		}
		out.close();
	}
}

#endif // MINISAT_LOGGER_H
