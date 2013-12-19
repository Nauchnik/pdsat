// class for solving SAT-problem using MPI

#ifndef mpi_solver_h
#define mpi_solver_h

#include "mpi_base.h"

class MPI_Solver : public MPI_Base
{
private:
	double *solving_times;
	vector<double> total_solving_times;
	unsigned solving_iteration_count;
	unsigned interrupted_count;
	string base_solving_info_file_name;
	string solving_info_file_name;
	double max_solving_time_koef;
	double finding_first_sat_time;
public:
	// Constructor/Destructor:
    MPI_Solver( );
    ~MPI_Solver( );
	
	int orig_tasks_count;
	int full_mask_tasks_count;
	int exch_activ;
	unsigned skip_tasks;
	double prev_med_time_sum;
	
	bool MPI_Solve( int argc, char **argv );
	bool ControlProcessSolve( );
	bool ComputeProcessSolve( );

	bool MPI_ConseqSolve( int argc, char **argv );
	bool WriteTimeToFile( double whole_time_sec );
	bool cpuTimeInHours( double full_seconds, int &real_hours, int &real_minutes, int &real_seconds );

	void CollectAssumptionsFiles();
	void PrintParams( );
	void WriteSolvingTimeInfo( double *solving_times, unsigned solved_tasks_count );
	void AddSolvingTimeToArray( ProblemStates cur_problem_state, double cnf_time_from_node, 
		                        double *solving_times );
	bool SolverRun( Solver *&S, unsigned &local_interrupted_count, int &process_sat_count, double &cnf_time_from_node, 
				    int current_task_index );
};

#endif
