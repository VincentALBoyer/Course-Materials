#pragma once

#include "../Solver.hpp"

class LPSolver :  public Solver {

protected:
	Env _env;

	bool _isIP;

	Solver::SolverStatus _status;

public:
	LPSolver(Instance* I) : Solver(I, "LP"), _env(I), _isIP(false), _status(Solver::SolverStatus::Unknown) {}
	~LPSolver() { }

	Env* env() { return &_env; }
	
	double gap() { return _env.getGap(); }
	Solution* recoversolution() {
		//if (_status != Solver::SolverStatus::Optimal && _status != Solver::SolverStatus::Feasible ) return NULL;
		return _env.recoverSolution(); 
	}

	NumVar* getVar(int k) { return _env.getVar(k); }
	virtual NumVar* getVar(FlowEdge* e);

	Solver::SolverStatus getStatus() { return _status; }

	void convertToIP();

	void convertToLP();

	void enableOutput(bool enable = true) {
		if (enable) _env.getCplex().setOut(std::cout);
		else _env.getCplex().setOut(_env.getNullStream());
	}
};

class LPCplexSolver : public LPSolver {

	void solvemethod(Solution* S) {
		_env.setTimLim(_timlim);
		_env.setThreads(_threads);
		_env.setMinGap(_mingap);
		_status = _env.runsolver();
	}

public:
	LPCplexSolver(Instance* I);
	~LPCplexSolver() {}

	
};

class CuttingPlaneSolver : public LPCplexSolver {

	bool _runonce;

	int _hasrun;

	std::vector<NumVar*> getsortedVars();

	std::vector<Constraint*> separateGomoryCuts(size_t SizeMax = std::numeric_limits<size_t>::max());

	std::vector<Constraint*> separateSECCuts(size_t SizeMax = std::numeric_limits<size_t>::max());


	void solvemethod(Solution* S) override;

public:
	CuttingPlaneSolver(Instance* I, bool runonce = true) : LPCplexSolver(I), _runonce(runonce), _hasrun(false) {
	
	
	}

	~CuttingPlaneSolver() {}

};



