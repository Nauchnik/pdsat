// class for solving SAT-problem using MPI

#ifndef mpi_solver_h
#define mpi_solver_h

#include "mpi_base.h"

struct satisfying_assignment
{
	std::string str;
	double solving_time;
};

class MPI_Solver : public MPI_Base
{
private:
	double solving_times[SOLVING_TIME_LEN];
	std::vector<double> total_solving_times;
	unsigned interrupted_count;
	std::string base_solving_info_file_name;
	double finding_first_sat_time;
	double total_start_time;
public:
	// Constructor/Destructor:
    MPI_Solver( );
    ~MPI_Solver( );
	
	unsigned solving_iteration_count;
	std::string solving_info_file_name;
	int full_mask_tasks_count;
	unsigned skip_tasks;
	double prev_med_time_sum;
	double max_solving_time_koef;
	bool no_increm;
	bool isCollectInterruptedInstances;
	int variables_each_integer;
	
	bool MPI_Solve( int argc, char **argv );
	bool ControlProcessSolve( std::vector<int> extern_var_choose_order, 
		                      std::vector<std::vector<bool>> &interrupted_problems_var_values,
							  std::vector<satisfying_assignment> &satisfying_assignments );
	bool ComputeProcessSolve();	
	
	bool MPI_ConseqSolve( int argc, char **argv );
	bool WriteTimeToFile( double whole_time_sec );

	void PrintParams();
	void WriteSolvingTimeInfo( double *solving_times, unsigned solved_tasks_count );
	void AddSolvingTimeToArray( ProblemStates cur_problem_state, double cnf_time_from_node, 
   	                            double *solving_times );
	bool SolverRun( Minisat::Solver *&S, unsigned long long &process_sat_count, double &cnf_time_from_node, 
				    int current_task_index, std::vector< std::vector<bool> > &interrupted_problems_var_values_from_process,
					std::vector< std::vector<bool> > &sat_assignment_from_process );
};

#endif
