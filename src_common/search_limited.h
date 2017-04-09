
//search limited
lbool Solver::search_limited() {

	assert(ok);
	int         backtrack_level;
	int         conflictC = 0;
	vec<Lit>    learnt_clause;
	starts++;
	//lbool status = l_Undef;
	// BOINC mode - added to speedup solving Latin square problems and decreasie using RAM
	
	for (;;) {		
		CRef confl = propagate();
		if (confl != CRef_Undef) {
			// CONFLICT
			conflicts++; conflictC++;
			if (decisionLevel() == 0) return l_False;

			learnt_clause.clear();
			analyze(confl, learnt_clause, backtrack_level);
			cancelUntil(backtrack_level);
			
			if (learnt_clause.size() == 1) {
				uncheckedEnqueue(learnt_clause[0]);
			}
			else {		
				CRef cr = ca.alloc(learnt_clause, true);
				learnts.push(cr);
				attachClause(cr);
				//claBumpActivity(ca[cr]);
				ca[cr].activity() = LBD(ca[cr]);
				uncheckedEnqueue(learnt_clause[0], cr);
			}

			varDecayActivity();
			claDecayActivity();

			if (--learntsize_adjust_cnt == 0) {
				learntsize_adjust_confl *= learntsize_adjust_inc;
				learntsize_adjust_cnt = (int)learntsize_adjust_confl;
				// max_learnts             *= learntsize_inc;
			}

		}
		else {
			// NO CONFLICT		
			/*
			// Simplify the set of problem clauses:
			if (decisionLevel() == 0 && !simplify())
				return l_False;

			if (learnts.size() - nAssigns() >= max_learnts)
				// Reduce the set of learnt clauses:
				reduceDB();
				*/
			Lit next = lit_Undef;
			while (decisionLevel() < assumptions.size()) {
				// Perform user provided assumption:
				Lit p = assumptions[decisionLevel()];
				if (value(p) == l_True) {
					// Dummy decision level:
					newDecisionLevel();
				}
				else if (value(p) == l_False) {					
					analyzeFinal(~p, conflict);					
					return l_False;
				}
				else {
					next = p;
					break;
				}
			}

			if (next == lit_Undef) {
				// New variable decision:
				return l_Undef;
			}

			// Increase decision level and enqueue 'next'
			newDecisionLevel();
			
			uncheckedEnqueue(next);
		}
	}
}