#include "minisat22_wrapper.h"
#include <iostream>

void minisat22_wrapper :: readClause(StreamBuffer& in, vec<Lit>& lits) 
{
    int     parsed_lit, var;
    lits.clear();
    for (;;){
        parsed_lit = parseInt(in);
        if (parsed_lit == 0) break;
        var = abs(parsed_lit)-1;
        lits.push( (parsed_lit > 0) ? mkLit(var) : ~mkLit(var) );
    }
}

void minisat22_wrapper :: parse_DIMACS_to_problem(std::istream& input, Problem& cnf)
{
	StreamBuffer in(input);
	Disjunct* lits = 0;
    int vars    = 0;
    int clauses = 0;
    int cnt     = 0;
    for (;;){
        skipWhitespace(in);
		if (in.eof()) break;
        else if (*in == 'p'){
            if (eagerMatch(in, "p cnf")){
                vars    = parseInt(in);
                clauses = parseInt(in);
            }else{
                printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
            }
        } else if (*in == 'c' || *in == 'p')
            skipLine(in);
        else{
            cnt++;
			lits = new Disjunct;
            readClause(in, *lits);
			cnf.push_back(lits);
		}
    }
}

// read cnf from resource file
void minisat22_wrapper :: parse_DIMACS_from_inc( std::vector<int> &cnf_vec, Problem &cnf )
{
	Disjunct *lits = 0;
	int parsed_lit, var;
	bool IsNewClause = true;
	
	for ( unsigned i = 0; i < cnf_vec.size(); i++ ) {
		parsed_lit = cnf_vec[i];
		if ( parsed_lit == 0 ) {
			cnf.push_back(lits);
			IsNewClause = true;
			continue;
		}
		var = abs(parsed_lit)-1;
		if ( IsNewClause ) {
			lits = new Disjunct;
			IsNewClause = false;
		}
		lits->push( (parsed_lit > 0) ? mkLit(var) : ~mkLit(var) );
    }
}

void minisat22_wrapper :: printProblem(const Problem& p, std::ostream& out)
{
	for(size_t i = 0; i < p.size(); i++) {
		for(int j = 0; j < p[i]->size(); j++) {
			Lit& lit = (*p[i])[j];
			out << (sign(lit)?"-":"") << var(lit)+1 << (j+1==p[i]->size()?" 0\n":" ");
		}
	}
}