
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

namespace Minisat
{

// Дескриптор для дапма состояния решателя
struct SolverStateDesc
{
	SolverStateDesc() {}

	// Mode of operation:

	int       verbosity;
	double    var_decay;
	double    clause_decay;
	double    random_var_freq;
	double    random_seed;
	bool      luby_restart;
	int       ccmin_mode;
	int       phase_saving;
	bool      rnd_pol;
	bool      rnd_init_act;
	double    garbage_frac;
	int       min_learnts_lim;

	int       restart_first;
	double    restart_inc;
	double    learntsize_factor;
	double    learntsize_inc;

	int       learntsize_adjust_start_confl;
	double    learntsize_adjust_inc;
	double    max_learnts;
	double    learntsize_adjust_confl;
	int       learntsize_adjust_cnt;

	// Solver state:

	bool      ok;
	double    cla_inc;
	double    var_inc;
	int       qhead;
	int       simpDB_assigns;
	int64_t   simpDB_props;
	double    progress_estimate;
	bool      remove_satisfied;
	Var       next_var;

	// Statistics:

	uint64_t  solves;
	uint64_t  starts;
	uint64_t  decisions;
	uint64_t  rnd_decisions;
	uint64_t  propagations;
	uint64_t  conflicts;
	uint64_t  dec_vars;
	uint64_t  num_clauses;
	uint64_t  num_learnts;
	uint64_t  clauses_literals;
	uint64_t  learnts_literals;
	uint64_t  max_literals;
	uint64_t  tot_literals;

	// Resource contraints:

	int64_t   conflict_budget;
	int64_t   propagation_budget;
	bool      asynch_interrupt;

	// Data collections:

	uint32_t  clauses_size;
	uint32_t  learnts_size;
	uint32_t  trail_size;
	uint32_t  trail_lim_size;
	uint32_t  assumptions_size;
	uint32_t  activity_size;
	uint32_t  assigns_size;
	uint32_t  polarity_size;
	uint32_t  user_pol_size;
	uint32_t  decision_size;
	uint32_t  vardata_size;
	uint32_t  watches_occs_size;
	uint32_t  watches_dirty_size;
	uint32_t  watches_dirties_size;
	uint32_t  order_heap_size;
	uint32_t  order_heap_indices_size;
	uint32_t  released_vars_size;
	uint32_t  free_vars_size;

	// Memory pool: 

	uint32_t  allocator_capacity;
	uint32_t  allocator_size;
	uint32_t  allocator_wasted;
};

class Solver;

class SolverStateAccessor
{
public:

	SolverStateAccessor(Solver& solver)
		: solver_(solver)
	{}

	void GetSolverStateDesc(SolverStateDesc& desc) const;

	void SetSolverStateDesc(const SolverStateDesc& desc);

	//! Сохранить состояние решателя в блоб
	void WriteStateBlob(const std::string& filename) const;

	//! Восстановить состояние решателя из блоба
	void ReadStateBlob(const std::string& filename);

	//! Сохранить состояние решателя в текстовом файле
	void WriteStateText(const std::string& filename) const;	

private:

	Solver& solver_;
};

std::ostream& operator<<(std::ostream& out, const Lit& lit);
std::ostream& operator<<(std::ostream& out, const lbool& lb);
std::ostream& operator<<(std::ostream& out, const SolverStateDesc& desc);

} // namespace Minisat
