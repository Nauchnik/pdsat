#include "../core/Solver.h"
#include "SolverStateAccessor.h"

namespace Minisat
{

namespace
{

template <typename T>
void WriteOut(std::ostream& out, std::string msg, const vec<T>& data)
{
	if(!msg.empty())
		out << std::endl << msg << std::endl;

	for(int i = 0; i < data.size(); ++i)
		out << data[i] << " ";
	out << std::endl;
}

template<class K, class V, class MkIndex>
void WriteOut(std::ostream& out, std::string msg, const IntMap<K, V, MkIndex>& data)
{
	if(!msg.empty())
		out << std::endl << msg << std::endl;

	for(const V* it = data.begin(); it != data.end(); ++it)
		out << *it << " ";
	out << std::endl;
}

template <class T>
void WriteBlob(std::ostream& out, const vec<T>& data)
{
	const std::size_t size = data.size() * sizeof(T);
	if(!size) return;
	const char* ptr = reinterpret_cast<const char*>((const T*)data);
	out.write(ptr, size);
}

template<class K, class V, class MkIndex>
void WriteBlob(std::ostream& out, const IntMap<K, V, MkIndex>& data)
{
	const std::size_t size = data.size() * sizeof(V);
	if(!size) return;
	const char* ptr = reinterpret_cast<const char*>(data.begin());
	out.write(ptr, size);
}

template <class T>
void ReadBlob(std::istream& in, std::size_t size, vec<T>& data)
{
	if(!size) return;
	data.growTo(size);
	const std::size_t read_size = size * sizeof(T);
	char* ptr = reinterpret_cast<char*>((T*)data);
	in.read(ptr, read_size);
}

template<class K, class V, class MkIndex>
void ReadBlob(std::istream& in, std::size_t size, IntMap<K, V, MkIndex>& data)
{
	if(!size) return;
	data.reserve(size);
	const std::size_t read_size = size * sizeof(V);
	char* ptr = reinterpret_cast<char*>(data.begin());
	in.read(ptr, read_size);
}

} // local namespace

void SolverStateAccessor::ReadStateBlob(const std::string& filename)
{
	std::ifstream in(filename.c_str(), std::ios::in | std::ios::binary);
	if(!in.is_open()) {
		//throw std::runtime_error("ReadStateBlob: can't open file " + filename);
		std::cerr << "ReadStateBlob: can't open file " << filename << std::endl;
		exit(1);
	}

	SolverStateDesc desc;
	in.read(reinterpret_cast<char*>(&desc), sizeof(desc));
	SetSolverStateDesc(desc);

	ReadBlob(in, desc.clauses_size, solver_.clauses);
	ReadBlob(in, desc.learnts_size, solver_.learnts);
	ReadBlob(in, desc.trail_size, solver_.trail);
	ReadBlob(in, desc.trail_lim_size, solver_.trail_lim);
	ReadBlob(in, desc.assumptions_size, solver_.assumptions);
	ReadBlob(in, desc.activity_size, solver_.activity);
	ReadBlob(in, desc.assigns_size, solver_.assigns);
	ReadBlob(in, desc.polarity_size, solver_.polarity);
	ReadBlob(in, desc.user_pol_size, solver_.user_pol);
	ReadBlob(in, desc.decision_size, solver_.decision);
	ReadBlob(in, desc.vardata_size, solver_.vardata);
	
	//! watches
	solver_.watches.occs.reserve(desc.watches_occs_size);
	for(vec<Solver::Watcher>* it = solver_.watches.occs.begin(); it != solver_.watches.occs.end(); ++it)
	{
		std::size_t size(0);
		in.read(reinterpret_cast<char*>(&size), sizeof(size));
		if(size)
		{
			ReadBlob(in, size, *it);
		}
	}
	ReadBlob(in, desc.watches_dirty_size, solver_.watches.dirty);
	ReadBlob(in, desc.watches_dirties_size, solver_.watches.dirties);

	//! heap
	ReadBlob(in, desc.order_heap_size, solver_.order_heap.heap);
	ReadBlob(in, desc.order_heap_indices_size, solver_.order_heap.indices);

	ReadBlob(in, desc.released_vars_size, solver_.released_vars);
	ReadBlob(in, desc.free_vars_size, solver_.free_vars);
	ReadBlob(in, desc.seen_size, solver_.seen);
	ReadBlob(in, desc.analyze_stack_size, solver_.analyze_stack);
	ReadBlob(in, desc.analyze_toclear_size, solver_.analyze_toclear);
	ReadBlob(in, desc.add_tmp_size, solver_.add_tmp);

	//! memory allocator
	{
		solver_.ca.ra.capacity(desc.allocator_capacity);
		const std::size_t memory_size = desc.allocator_capacity * ClauseAllocator::Unit_Size;
		char* memory_ptr = reinterpret_cast<char*>(solver_.ca.ra.memory);
		in.read(memory_ptr, memory_size);
		//assert(solver_.ca.ra.cap == desc.allocator_capacity);
		solver_.ca.ra.sz = desc.allocator_size;
		solver_.ca.ra.wasted_ = desc.allocator_wasted;
		solver_.ca.ra.cap = desc.allocator_capacity;
	}
}

void SolverStateAccessor::WriteStateBlob(const std::string& filename) const
{
	std::ofstream out(filename.c_str(), std::ios::out | std::ios::binary);
	if(!out.is_open()) {
		//throw std::runtime_error("WriteStateBlob: can't open file " + filename);
		std::cerr << "ReadStateBlob: can't open file " << filename << std::endl;
		exit(1);
	}

	SolverStateDesc desc;
	GetSolverStateDesc(desc);
	out.write(reinterpret_cast<char*>(&desc), sizeof(desc));

	WriteBlob(out, solver_.clauses);
	WriteBlob(out, solver_.learnts);
	WriteBlob(out, solver_.trail);
	WriteBlob(out, solver_.trail_lim);
	WriteBlob(out, solver_.assumptions);
	WriteBlob(out, solver_.activity);
	WriteBlob(out, solver_.assigns);
	WriteBlob(out, solver_.polarity);
	WriteBlob(out, solver_.user_pol);
	WriteBlob(out, solver_.decision);
	WriteBlob(out, solver_.vardata);
	
	//! watches
	for(vec<Solver::Watcher>* it = solver_.watches.occs.begin(); it != solver_.watches.occs.end(); ++it)
	{
		const std::size_t vec_size = it->size();
		out.write(reinterpret_cast<const char*>(&vec_size), sizeof(vec_size));
		WriteBlob(out, *it);
	}
	WriteBlob(out, solver_.watches.dirty);
	WriteBlob(out, solver_.watches.dirties);
	
	//! heap
	WriteBlob(out, solver_.order_heap.heap);
	WriteBlob(out, solver_.order_heap.indices);

	WriteBlob(out, solver_.released_vars);
	WriteBlob(out, solver_.free_vars);
	WriteBlob(out, solver_.seen);
	WriteBlob(out, solver_.analyze_stack);
	WriteBlob(out, solver_.analyze_toclear);
	WriteBlob(out, solver_.add_tmp);

	//! memory allocator
	{
		//! Сохряняем всю выделенную память
		const std::size_t memory_size = solver_.ca.ra.capacity() * ClauseAllocator::Unit_Size;
		const char* memory_ptr = reinterpret_cast<const char*>(solver_.ca.ra.memory);
		out.write(memory_ptr, memory_size);
	}
}

void SolverStateAccessor::GetSolverStateDesc(SolverStateDesc& desc) const
{
	desc.verbosity = solver_.verbosity;
	desc.var_decay = solver_.var_decay;
	desc.clause_decay = solver_.clause_decay;
	desc.random_var_freq = solver_.random_var_freq;
	desc.random_seed = solver_.random_seed;
	desc.luby_restart = solver_.luby_restart;
	desc.ccmin_mode = solver_.ccmin_mode;
	desc.phase_saving = solver_.phase_saving;
	desc.rnd_pol = solver_.rnd_pol;
	desc.rnd_init_act = solver_.rnd_init_act;
	desc.garbage_frac = solver_.garbage_frac;
	desc.min_learnts_lim = solver_.min_learnts_lim;
	desc.restart_first = solver_.restart_first;
	desc.restart_inc = solver_.restart_inc;
	desc.learntsize_factor = solver_.learntsize_factor;
	desc.learntsize_inc = solver_.learntsize_inc;
	desc.learntsize_adjust_start_confl = solver_.learntsize_adjust_start_confl;
	desc.learntsize_adjust_inc = solver_.learntsize_adjust_inc;
	desc.max_learnts = solver_.max_learnts;
	desc.learntsize_adjust_confl = solver_.learntsize_adjust_confl;
	desc.learntsize_adjust_cnt = solver_.learntsize_adjust_cnt;
	desc.ok = solver_.ok;
	desc.cla_inc = solver_.cla_inc;
	desc.var_inc = solver_.var_inc;
	desc.qhead = solver_.qhead;
	desc.simpDB_assigns = solver_.simpDB_assigns;
	desc.simpDB_props = solver_.simpDB_props;
	desc.progress_estimate = solver_.progress_estimate;
	desc.remove_satisfied = solver_.remove_satisfied;
	desc.next_var = solver_.next_var;
	desc.solves = solver_.solves;
	desc.starts = solver_.starts;
	desc.decisions = solver_.decisions;
	desc.rnd_decisions = solver_.rnd_decisions;
	desc.propagations = solver_.propagations;
	desc.conflicts = solver_.conflicts;
	desc.dec_vars = solver_.dec_vars;
	desc.num_clauses = solver_.num_clauses;
	desc.num_learnts = solver_.num_learnts;
	desc.clauses_literals = solver_.clauses_literals;
	desc.learnts_literals = solver_.learnts_literals;
	desc.max_literals = solver_.max_literals;
	desc.tot_literals = solver_.tot_literals;
	desc.conflict_budget = solver_.conflict_budget;
	desc.propagation_budget = solver_.propagation_budget;
	desc.asynch_interrupt = solver_.asynch_interrupt;
	desc.clauses_size = solver_.clauses.size();
	desc.learnts_size = solver_.learnts.size();
	desc.trail_size = solver_.trail.size();
	desc.trail_lim_size = solver_.trail_lim.size();
	desc.assumptions_size = solver_.assumptions.size();
	desc.activity_size = solver_.activity.size();
	desc.assigns_size = solver_.assigns.size();
	desc.polarity_size = solver_.polarity.size();
	desc.user_pol_size = solver_.user_pol.size();
	desc.decision_size = solver_.decision.size();
	desc.vardata_size = solver_.vardata.size();
	desc.watches_occs_size = solver_.watches.occs.size();
	desc.watches_dirty_size = solver_.watches.dirty.size();
	desc.watches_dirties_size = solver_.watches.dirties.size();
	desc.order_heap_size = solver_.order_heap.size();
	desc.order_heap_indices_size = solver_.order_heap.indices.size();
	desc.released_vars_size = solver_.released_vars.size();
	desc.free_vars_size = solver_.free_vars.size();
	desc.seen_size = solver_.seen.size();
	desc.analyze_stack_size = solver_.analyze_stack.size();
	desc.analyze_toclear_size = solver_.analyze_toclear.size();
	desc.add_tmp_size = solver_.add_tmp.size();
	desc.allocator_capacity = solver_.ca.ra.capacity();
	desc.allocator_size = solver_.ca.ra.size();
	desc.allocator_wasted = solver_.ca.ra.wasted();
}

void SolverStateAccessor::SetSolverStateDesc(const SolverStateDesc& desc)
{
	solver_.verbosity = desc.verbosity;
	solver_.var_decay = desc.var_decay;
	solver_.clause_decay = desc.clause_decay;
	solver_.random_var_freq = desc.random_var_freq;
	solver_.random_seed = desc.random_seed;
	solver_.luby_restart = desc.luby_restart;
	solver_.ccmin_mode = desc.ccmin_mode;
	solver_.phase_saving = desc.phase_saving;
	solver_.rnd_pol = desc.rnd_pol;
	solver_.rnd_init_act = desc.rnd_init_act;
	solver_.garbage_frac = desc.garbage_frac;
	solver_.min_learnts_lim = desc.min_learnts_lim;
	solver_.restart_first = desc.restart_first;
	solver_.restart_inc = desc.restart_inc;
	solver_.learntsize_factor = desc.learntsize_factor;
	solver_.learntsize_inc = desc.learntsize_inc;
	solver_.learntsize_adjust_start_confl = desc.learntsize_adjust_start_confl;
	solver_.learntsize_adjust_inc = desc.learntsize_adjust_inc;
	solver_.max_learnts = desc.max_learnts;
	solver_.learntsize_adjust_confl = desc.learntsize_adjust_confl;
	solver_.learntsize_adjust_cnt = desc.learntsize_adjust_cnt;
	solver_.ok = desc.ok;
	solver_.cla_inc = desc.cla_inc;
	solver_.var_inc = desc.var_inc;
	solver_.qhead = desc.qhead;
	solver_.simpDB_assigns = desc.simpDB_assigns;
	solver_.simpDB_props = desc.simpDB_props;
	solver_.progress_estimate = desc.progress_estimate;
	solver_.remove_satisfied = desc.remove_satisfied;
	solver_.next_var = desc.next_var;
	solver_.solves = desc.solves;
	solver_.starts = desc.starts;
	solver_.decisions = desc.decisions;
	solver_.rnd_decisions = desc.rnd_decisions;
	solver_.propagations = desc.propagations;
	solver_.conflicts = desc.conflicts;
	solver_.dec_vars = desc.dec_vars;
	solver_.num_clauses = desc.num_clauses;
	solver_.num_learnts = desc.num_learnts;
	solver_.clauses_literals = desc.clauses_literals;
	solver_.learnts_literals = desc.learnts_literals;
	solver_.max_literals = desc.max_literals;
	solver_.tot_literals = desc.tot_literals;
	solver_.conflict_budget = desc.conflict_budget;
	solver_.propagation_budget = desc.propagation_budget;
	solver_.asynch_interrupt = desc.asynch_interrupt;
}

void SolverStateAccessor::WriteStateText(const std::string& filename) const
{
	std::ofstream out(filename.c_str(), std::ios::out);
	if(!out.is_open()) {
		//throw std::runtime_error("WriteStateText: can't open file " + filename);
		std::cerr << "ReadStateBlob: can't open file " << filename << std::endl;
		exit(1);
	}
	
	SolverStateDesc desc;
	GetSolverStateDesc(desc);
	out << desc << std::endl;

	WriteOut(out, "clauses:", solver_.clauses);
	WriteOut(out, "learnts:", solver_.learnts);
	WriteOut(out, "trail:", solver_.trail);
	WriteOut(out, "trail_lim:", solver_.trail_lim);
	WriteOut(out, "assumptions:", solver_.assumptions);
	WriteOut(out, "activity:", solver_.activity);
	WriteOut(out, "assigns:", solver_.assigns);
	WriteOut(out, "polarity:", solver_.polarity);
	WriteOut(out, "user_pol:", solver_.user_pol);
	WriteOut(out, "decision:", solver_.decision);
	WriteOut(out, "vardata:", solver_.vardata);
	out << std::endl << "watches.occs:" << std::endl;
	for(vec<Solver::Watcher>* it = solver_.watches.occs.begin(); it != solver_.watches.occs.end(); ++it)
	{
		WriteOut(out, "", *it);
	}
	WriteOut(out, "watches.dirty:", solver_.watches.dirty);
	WriteOut(out, "watches.dirties:", solver_.watches.dirties);
	WriteOut(out, "order_heap.heap:", solver_.order_heap.heap);
	WriteOut(out, "order_heap.indices:", solver_.order_heap.indices);
	WriteOut(out, "released_vars:", solver_.released_vars);
	WriteOut(out, "free_vars:", solver_.free_vars);
	WriteOut(out, "seen:", solver_.seen);
	WriteOut(out, "analyze_stack:", solver_.analyze_stack);
	WriteOut(out, "analyze_toclear:", solver_.analyze_toclear);
	WriteOut(out, "add_tmp:", solver_.add_tmp);
}

std::ostream& operator<<(std::ostream& out, const SolverStateDesc& desc)
{
	out << "verbosity: " << desc.verbosity << std::endl;
	out << "var_decay: " << desc.var_decay << std::endl;
	out << "clause_decay: " << desc.clause_decay << std::endl;
	out << "random_var_freq: " << desc.random_var_freq << std::endl;
	out << "random_seed: " << desc.random_seed << std::endl;
	out << "luby_restart: " << desc.luby_restart << std::endl;
	out << "ccmin_mode: " << desc.ccmin_mode << std::endl;
	out << "phase_saving: " << desc.phase_saving << std::endl;
	out << "rnd_pol: " << desc.rnd_pol << std::endl;
	out << "rnd_init_act: " << desc.rnd_init_act << std::endl;
	out << "garbage_frac: " << desc.garbage_frac << std::endl;
	out << "min_learnts_lim: " << desc.min_learnts_lim << std::endl;
	out << "restart_first: " << desc.restart_first << std::endl;
	out << "restart_inc: " << desc.restart_inc << std::endl;
	out << "learntsize_factor: " << desc.learntsize_factor << std::endl;
	out << "learntsize_inc: " << desc.learntsize_inc << std::endl;
	out << "learntsize_adjust_start_confl: " << desc.learntsize_adjust_start_confl << std::endl;
	out << "learntsize_adjust_inc: " << desc.learntsize_adjust_inc << std::endl;
	out << "max_learnts: " << desc.max_learnts << std::endl;
	out << "learntsize_adjust_confl: " << desc.learntsize_adjust_confl << std::endl;
	out << "learntsize_adjust_cnt: " << desc.learntsize_adjust_cnt << std::endl;
	out << "ok: " << desc.ok << std::endl;
	out << "cla_inc: " << desc.cla_inc << std::endl;
	out << "var_inc: " << desc.var_inc << std::endl;
	out << "qhead: " << desc.qhead << std::endl;
	out << "simpDB_assigns: " << desc.simpDB_assigns << std::endl;
	out << "simpDB_props: " << desc.simpDB_props << std::endl;
	out << "progress_estimate: " << desc.progress_estimate << std::endl;
	out << "remove_satisfied: " << desc.remove_satisfied << std::endl;
	out << "next_var: " << desc.next_var << std::endl;
	out << "solves: " << desc.solves << std::endl;
	out << "starts: " << desc.starts << std::endl;
	out << "decisions: " << desc.decisions << std::endl;
	out << "rnd_decisions: " << desc.rnd_decisions << std::endl;
	out << "propagations: " << desc.propagations << std::endl;
	out << "conflicts: " << desc.conflicts << std::endl;
	out << "dec_vars: " << desc.dec_vars << std::endl;
	out << "num_clauses: " << desc.num_clauses << std::endl;
	out << "num_learnts: " << desc.num_learnts << std::endl;
	out << "clauses_literals: " << desc.clauses_literals << std::endl;
	out << "learnts_literals: " << desc.learnts_literals << std::endl;
	out << "max_literals: " << desc.max_literals << std::endl;
	out << "tot_literals: " << desc.tot_literals << std::endl;
	out << "conflict_budget: " << desc.conflict_budget << std::endl;
	out << "propagation_budget: " << desc.propagation_budget << std::endl;
	out << "asynch_interrupt: " << desc.asynch_interrupt << std::endl;
	out << "clauses_size: " << desc.clauses_size << std::endl;
	out << "learnts_size: " << desc.learnts_size << std::endl;
	out << "trail_size: " << desc.trail_size << std::endl;
	out << "trail_lim_size: " << desc.trail_lim_size << std::endl;
	out << "assumptions_size: " << desc.assumptions_size << std::endl;
	out << "activity_size: " << desc.activity_size << std::endl;
	out << "assigns_size: " << desc.assigns_size << std::endl;
	out << "polarity_size: " << desc.polarity_size << std::endl;
	out << "user_pol_size: " << desc.user_pol_size << std::endl;
	out << "decision_size: " << desc.decision_size << std::endl;
	out << "vardata_size: " << desc.vardata_size << std::endl;
	out << "watches_occs_size: " << desc.watches_occs_size << std::endl;
	out << "watches_dirty_size: " << desc.watches_dirty_size << std::endl;
	out << "watches_dirties_size: " << desc.watches_dirties_size << std::endl;
	out << "order_heap_size: " << desc.order_heap_size << std::endl;
	out << "order_heap_indices_size: " << desc.order_heap_indices_size << std::endl;
	out << "released_vars_size: " << desc.released_vars_size << std::endl;
	out << "free_vars_size: " << desc.free_vars_size << std::endl;
	out << "seen_size: " << desc.seen_size << std::endl;
	out << "analyze_stack_size: " << desc.analyze_stack_size << std::endl;
	out << "analyze_toclear_size: " << desc.analyze_toclear_size << std::endl;
	out << "add_tmp_size: " << desc.add_tmp_size << std::endl;
	out << "allocator_capacity: " << desc.allocator_capacity << std::endl;
	out << "allocator_size: " << desc.allocator_size << std::endl;
	out << "allocator_wasted: " << desc.allocator_wasted << std::endl;
	return out;
}

std::ostream& operator<<(std::ostream& out, const Lit& lit)
{
	out << lit.x;
	return out;
}

std::ostream& operator<<(std::ostream& out, const lbool& lb)
{
	out << toInt(lb);
	return out;
}

std::ostream& operator<<(std::ostream& out, const Solver::VarData& data)
{
	out << "(" << data.reason << "," << data.level << ")";
	return out;
}

std::ostream& operator<<(std::ostream& out, const Solver::Watcher& watcher)
{
	out << "(" << watcher.cref << "," << watcher.blocker << ")";
	return out;
}

std::ostream& operator<<(std::ostream& out, const Solver::ShrinkStackElem& elem)
{
	out << "(" << elem.i << "," << elem.l << ")";
	return out;
}

} // namespace Minisat
