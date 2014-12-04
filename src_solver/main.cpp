#include <errno.h>

#include <signal.h>
//#include <zlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "./minisat/utils/System.h"
#include "./minisat/utils/ParseUtils.h"
#include "./minisat/utils/Options.h"
#include "./minisat/core/Dimacs.h"
#include "./minisat/core/Solver.h"

using namespace Minisat;

static Solver* solver;

//#define USE_PROBLEM

// to do: добавить реализацию функци printStats()

///////////////////////////////////////////////////////////////////////////////
// реализация возможности добавления в решатель уже прочитанной КНФ

//#ifdef USE_PROBLEM

static void readClause(StreamBuffer& in, vec<Lit>& lits) {
    int     parsed_lit, var;
    lits.clear();
    for (;;){
        parsed_lit = parseInt(in);
        if (parsed_lit == 0) break;
        var = abs(parsed_lit)-1;
        lits.push( (parsed_lit > 0) ? mkLit(var) : ~mkLit(var) );
    }
}
static void parse_DIMACS(std::istream& input, Problem& cnf)
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
static void printProblem(const Problem& p, std::ostream& out)
{
	for(size_t i = 0; i < p.size(); i++)
	{
		for(int j = 0; j < p[i]->size(); j++)
		{
			Lit& lit = (*p[i])[j];
			out << (sign(lit)?"-":"") << var(lit)+1 << (j+1==p[i]->size()?" 0\n":" ");
		}
	}
}

int main(int argc, char* argv[])
{
#ifdef _DEBUG
	argc = 3;
	argv[1] = "bivium_template.cnf";
	argv[2] = "out";
#endif
	double time = cpuTime();
	std::ifstream in(argv[1], std::ios::in);
	if (!in.is_open())
        printf("ERROR! Could not open file: %s\n", argc == 1 ? "<stdin>" : argv[1]), exit(1);

	double initial_time = cpuTime();

	// КНФ считываем 1 раз
	Problem cnf;
	parse_DIMACS(in, cnf);
	in.close();
    double parsed_time = cpuTime();

	// передача КНФ решателю, можно делать в цикле
	Solver S;    
	S.verbosity = 2;
	S.addProblem(cnf);
	
    FILE* res = (argc >= 3) ? fopen(argv[2], "wb") : NULL;

	if (S.verbosity > 0){
        printf("|  Number of variables:  %12d                                         |\n", S.nVars());
        printf("|  Number of clauses:    %12d                                         |\n", S.nClauses()); }
        
    if (S.verbosity > 0){
        printf("|  Parse time:           %12.2f s                                       |\n", parsed_time - initial_time);
        printf("|                                                                       |\n"); }

	if (!S.simplify()){
        if (res != NULL) fprintf(res, "UNSAT\n"), fclose(res);
        if (S.verbosity > 0){
            printf("===============================================================================\n");
            printf("Solved by unit propagation\n");
            //printStats(S);
            printf("\n"); }
        printf("UNSATISFIABLE\n");
		printf( "time %f", cpuTime() - time);
        exit(20);
    }

	printf("|  after simplify:           %12.2f s                                       |\n", parsed_time - initial_time);

	//S.start_activity = 1;
	//S.resetIntervalVarActivity( 578, 777 );
	S.print_learnts = true;
    
    vec<Lit> dummy;
    lbool ret = S.solveLimited(dummy);
    if (S.verbosity > 0){
        //printStats(S);
        printf("\n"); }
    //printf(ret == l_True ? "SATISFIABLE\n" : ret == l_False ? "UNSATISFIABLE\n" : "INDETERMINATE\n");

    if (res != NULL){
        if (ret == l_True){
            fprintf(res, "SAT\n");
            for (int i = 0; i < S.nVars(); i++)
                if (S.model[i] != l_Undef)
                    fprintf(res, "%s%s%d", (i==0)?"":" ", (S.model[i]==l_True)?"":"-", i+1);
            fprintf(res, " 0\n");
        }else if (ret == l_False)
            fprintf(res, "UNSAT\n");
        else
            fprintf(res, "INDET\n");
        fclose(res);
    }
	printf( "time %f", cpuTime() - time);
    
#ifdef NDEBUG
        exit(ret == l_True ? 10 : ret == l_False ? 20 : 0);     // (faster than "return", which will invoke the destructor for 'Solver')
#else
        return (ret == l_True ? 10 : ret == l_False ? 20 : 0);
#endif
	
	//debug
	//std::ofstream out("out.cnf", std::ios::out);
	//printProblem(cnf, out);
	//out.close();
	
	return 0;
}
/*
#else

///////////////////////////////////////////////////////////////////////////////
// стандартное использование - загрузка КНФ из файла

int main(int argc, char* argv[])
{
#ifdef _DEBUG
	argc = 2;
	argv[1] = "bivium_template.cnf";
#endif
	try
	{
		setUsageHelp("USAGE: %s [options] <input-file-cnf> <input-file-assigns> <result-output-file>\n\n  where input may be either in plain or gzipped DIMACS.\n");
		// Extra options:
        //
        IntOption    verb   ("MAIN", "verb",   "Verbosity level (0=silent, 1=some, 2=more).", 0, IntRange(0, 2));
        IntOption    cpu_lim("MAIN", "cpu-lim","Limit on CPU time allowed in seconds.\n", INT32_MAX, IntRange(0, INT32_MAX));
        IntOption    mem_lim("MAIN", "mem-lim","Limit on memory usage in megabytes.\n", INT32_MAX, IntRange(0, INT32_MAX));
        
        parseOptions(argc, argv, true);

        Solver S;
        double initial_time = cpuTime();

        S.verbosity = verb;
        
        solver = &S;

		if (argc == 1)
            printf("Reading from standard input... Use '--help' for help.\n");
        
        //gzFile in = (argc == 1) ? gzdopen(0, "rb") : gzopen(argv[1], "rb");
		std::ifstream in(argv[1], std::ios::in);
		if (!in.is_open())
            printf("ERROR! Could not open file: %s\n", argc == 1 ? "<stdin>" : argv[1]), exit(1);

		//S.verbosity = 1;

        if (S.verbosity > 0){
            printf("============================[ Problem Statistics ]=============================\n");
            printf("|                                                                             |\n"); }
 
    	//S.size_actitity = 0;
    	if(argc >= 5) // was 4
    	{
			// прочитать файл с информацией об активности переменных (?)
    		std::ifstream input_activ(argv[3], std::ios::in);
			if(!input_activ.is_open())
				printf("ERROR! Could not open activity file: %s\n", argv[3]), exit(1);
    		//read_Activity(input_activ, S);
    	}        
        
       // parse_DIMACS(in, S);

        //gzclose(in);
		in.close();

        FILE* res = (argc >= 3) ? fopen(argv[2], "wb") : NULL;

		if (S.verbosity > 0){
            printf("|  Number of variables:  %12d                                         |\n", S.nVars());
            printf("|  Number of clauses:    %12d                                         |\n", S.nClauses()); }
        
        double parsed_time = cpuTime();
        if (S.verbosity > 0){
            printf("|  Parse time:           %12.2f s                                       |\n", parsed_time - initial_time);
            printf("|                                                                             |\n"); }

		if (!S.simplify()){
            if (res != NULL) fprintf(res, "UNSAT\n"), fclose(res);
            if (S.verbosity > 0){
                printf("===============================================================================\n");
                printf("Solved by unit propagation\n");
                //printStats(S);
                printf("\n"); }
            printf("UNSATISFIABLE\n");
            exit(20);
        }
 
		vec<Lit> dummy;

        lbool ret = S.solveLimited(dummy);
        if (S.verbosity > 0){
            printf("\n"); }
        printf(ret == l_True ? "SATISFIABLE\n" : ret == l_False ? "UNSATISFIABLE\n" : "INDETERMINATE\n");
        if (res != NULL){
            if (ret == l_True){
                fprintf(res, "SAT\n");
                for (int i = 0; i < S.nVars(); i++)
                    if (S.model[i] != l_Undef)
                        fprintf(res, "%s%s%d", (i==0)?"":" ", (S.model[i]==l_True)?"":"-", i+1);
                fprintf(res, " 0\n");
            }else if (ret == l_False)
                fprintf(res, "UNSAT\n");
            else
                fprintf(res, "INDET\n");
            fclose(res);
        }
		double solved_time = cpuTime();
        printf("|  Solve time:           %12.2f s                                       |\n", solved_time - initial_time);
        
#ifdef NDEBUG
        exit(ret == l_True ? 10 : ret == l_False ? 20 : 0);     // (faster than "return", which will invoke the destructor for 'Solver')
#else
        return (ret == l_True ? 10 : ret == l_False ? 20 : 0);
#endif
	}
	catch (OutOfMemoryException&){
        printf("===============================================================================\n");
        printf("INDETERMINATE\n");
        exit(0);
    }
	//Solver solver;
	return 0;
}

#endif*/