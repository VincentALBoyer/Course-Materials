#pragma once
#include "Solution.hpp"



typedef int param_t;

class Solver {

protected:
	Instance* _I;
	string _name;
	double _timlim;
	double _mingap;
	int _threads;
	chrono::time_point<chrono::system_clock> _start;
	chrono::time_point<chrono::system_clock> _end;

	void startTimer() { _start = chrono::system_clock::now(); } 

	void stopTimer() { _end = chrono::system_clock::now(); }

	virtual void solvemethod(Solution* S) = 0;

	string PadNumber(int num, int w = 2, char c = '0');
	string RedString(string s, int w = 20);

public:
	//Status of the solver
	enum class NodeStatus { PendingNode, InfeasibleNode, FeasibleNode, PrunedNode, IntegerNode, IntegerInfeasibleNode };
	enum class VarStatus { Down, Up, Free };
	enum SolverStatus { Optimal, Feasible, Infeasible, Unknown };
	enum ConstraintType { OnNode, OnEdge, Gomory, NotSpecified };
	enum VarType { EdgeVar, UndEdgeVar, SlackVar, CapVar };

	//parameters
	static const param_t TimLim = 0;
	static const param_t Gap = 1;
	static const param_t Threads = 2;

	Solver(Instance* I, string name);

	virtual ~Solver() {};

	Instance* instance() { return _I; }

	void setparam(param_t PARAM, double value);

	string name() { return _name; }

	double cpuTime() { return (double)chrono::duration_cast<std::chrono::seconds>(_end - _start).count(); }

	double timeLeft() { return _timlim - (double)chrono::duration_cast<std::chrono::seconds>(chrono::system_clock::now() - _start).count();	}
	
	void solve(Solution* S=NULL);

	virtual double gap() = 0;

	virtual Solution* recoversolution()=0;
    
    void save(string outputfile="None");

};

struct Edge {
	int i;
	int j;
	double w;

	bool operator==(const Edge& other) const { return i == other.i && j == other.j; }

	// Define the order relationship
	bool operator<(const Edge& other) const {
		if (i != other.i) return i < other.i;
		if (j != other.j) return j < other.j;
		return w < other.w;
	}
};

class NumVar;
class NumColumn;
class Constraint;

struct VarKey {
	int i, j;
	Solver::VarType type;
	bool operator==(const VarKey& other) const {
		return i == other.i && j == other.j && type == other.type;
	}
};

namespace std {
	template <>
	struct hash<VarKey> {
		std::size_t operator()(const VarKey& k) const {
			return ((std::hash<int>()(k.i) ^ (std::hash<int>()(k.j) << 1)) >> 1)
				^ (std::hash<int>()(static_cast<int>(k.type)) << 1);
		}
	};
}

class Env : public IloEnv {
protected: 
	friend class Solver;
	Instance* _I;

	IloModel _model;
	IloCplex _cplex;
	IloObjective _obj;

	vector<NumVar*> _X;
	vector<Constraint*> _constraints;
	std::unordered_map<VarKey, size_t> _varIndexMap;

	NumVar* getVar(int i, int j, Solver::VarType type) const;



public:

	//DELTA for Cut Acceptance
	const static double DELTA;

	Env(Instance* I) : IloEnv(), _I(I) {
		_model = IloModel(static_cast<IloEnv&>(*this));
		_cplex = IloCplex(_model);
		_obj = IloObjective(static_cast<IloEnv&>(*this), 0, IloObjective::Minimize);
		_model.add(_obj);
	}
	Env() : IloEnv(), _I(NULL) {}
	~Env();

	IloCplex& getCplex() { return _cplex; }

	void setInstance(Instance* I) { _I = I; }
	Instance* getInstance() const { return _I; }

	void addConstraint(Constraint* constraint);
	void addConstraints(vector<Constraint*>& constraints);
	void setupConstraint(Constraint* constraint);
	vector<Constraint*>& getConstraints() const { return const_cast<vector<Constraint*>&>(_constraints);}
	void removeConstraints(vector<Constraint*>& constraints);

	int nvars() { return (int)_X.size(); }
	void addVar(NumVar* X);
	//void addVars(vector<NumVar*>& X) { _X.insert(_X.end(), X.begin(), X.end()); }
    vector<NumVar*>& getVars() const { return const_cast<vector<NumVar*>&>(_X); }
	NumVar* getVar(int k) const { return _X[k]; }
	NumVar* getEdgeVar(int i, int j) const { return getVar(i, j, Solver::EdgeVar); }
	NumVar* getCapVar(int i) const { return getVar(i, -1, Solver::CapVar); }
	

	double getObjVal() const { return _cplex.getObjValue(); }
	double getVarVal(NumVar& X) const;
	IloNumArray getVarVals(std::vector<NumVar*>& X);
	IloNumArray getVarVals() { return getVarVals(_X); }

	IloNumArray getDuals();
	
	int getNumVars() const { return int(_X.size()); }
	int getNumConstraints() const { return int(_constraints.size()); }
	Constraint* getConstraint(int i) { return _constraints[i]; }

	//bool isSlackVar(int k);
	//bool isEdgeVar(int k);

	vector<double> getInvRow(int k);
	vector<double> getInvBasisRow(int k);
	double getInvRHS(int k);
	pair<int*, double*> getBHead();
	IloCplex::BasisStatusArray getBasisStatus(const std::vector<NumVar*>& X);
    std::vector<std::string> getSortedColNames();
	//We compute the reduced cost of a variable
	double getReducedCost(NumVar* x);
	

	void addCol(NumColumn* col);

	//We compute the distance to a cut
	double getDistanceToCut(Constraint& cut);

	//Criteria for Cut acceptance
	bool acceptCut(Constraint& cut) { return getDistanceToCut(cut) > Env::DELTA;}

	void setObjVarCoef(NumVar* X, double obj);
	void setObjective(std::function<double(NumVar*)> costfunction);
	void setTimLim(double timlim) { _cplex.setParam(IloCplex::TiLim, timlim); }
	void setMinGap(double gap) { _cplex.setParam(IloCplex::EpGap, gap); }
	void setThreads(int threads) { _cplex.setParam(IloCplex::Threads, threads); }
	double getGap();
	Solver::SolverStatus runsolver();

	Solution* recoverSolution();

	void convertToIP();
	void convertToLP();

	void exportModel(string filename) { _cplex.exportModel(filename.c_str()); }	

	void DebugCheckConstraints(const std::vector<Constraint*>& C);
	std::vector<double> extractVarVals(Solution* S) const;

};


class NumVar : public IloNumVar {

protected:
	static unsigned int _index;

	unsigned int _id;

	int _i;
	int _j;

	static string nameSlackGen(int id) { return "S_" + to_string(id); }

	static string nameEdgeGen(int i, int j);

	void consolidate();

	Solver::VarType _type;

	Env* _env;
	bool _ended;
public:



	//static unsigned int hash(int i, int j) { return i * 1000 + j; }

	NumVar(Env& env) : IloNumVar(static_cast<IloEnv&>(env), 0, IloInfinity), _id(_index++), _env(&env), _type(Solver::SlackVar), _i(_id), _j(_id), _ended(false) {
		consolidate();
	}
	NumVar(Env& env, int i, int j) : IloNumVar(static_cast<IloEnv&>(env), 0, 1), _id(_index++), _env(&env), _type(Solver::EdgeVar), _i(i), _j(j), _ended(false) {
		consolidate();
	}
	NumVar(Env& env, int i) : IloNumVar(static_cast<IloEnv&>(env), 0, IloInfinity), _id(_index++), _env(&env), _type(Solver::CapVar), _i(i), _j(-1), _ended(false) {
		consolidate();
	}
	NumVar(NumVar& var) : IloNumVar(var), _id(var._id), _env(var._env), _type(var._type), _i(var._i), _j(var._j), _ended(false) {}
	//NumVar(Env& env, IloNumColumn& col) : IloNumVar(col, 0, IloInfinity), _id(_index++), _env(&env), _type(Env::SlackVar) {
	//	consolidate();
	//}
	NumVar(Env& env, IloNumColumn& col, int i, int j) : IloNumVar(col, 0, 1), _id(_index++), _env(&env), _type(Solver::EdgeVar), _i(i), _j(j), _ended(false) {
		consolidate();
	}
	~NumVar() { 
		if(!_ended) this->end(); 
	}

	unsigned int getId() const { return _id; }

	Solver::VarType getType() const { return _type; }


	int geti() const { return _i; }
	int getj() const { return _j; }

	void setasended() { _ended = true; } // Set the variable as ended, so it won't be ended again in the destructor

};


class NumColumn : public std::map<int, double> {

	NumVar* _var;

	double _obj;
	int _i;
	int _j;

public:

	NumColumn(int i, int j, double obj) : std::map<int, double>(), _i(i), _j(j), _obj(obj), _var(NULL) {
		if (_i < _j) std::swap(_i, _j);
	}
	~NumColumn() {}

	double getObj() const { return _obj; }
	int geti() const { return _i; }
	int getj() const { return _j; }
	NumVar* getVar() const { return _var; }

	void setObj(double obj) { _obj = obj; }
	void add(int i, double val) {
		if (fabs(val) < EPSILON) return;
		if (isInModel()) return;

		this->insert_or_assign(i, val);
	}

	void setLB(double lb) { _var->setLB(lb); }
	void setUB(double ub) { _var->setUB(ub); }

	bool isInModel() const { return _var != NULL; }
	void setInModel(NumVar* x) {
		//this->clear(); //We do not need the information anymore (memory optimization)
		_var = x;
	}

	double getReducedCost(IloNumArray& duals) {
		return std::accumulate(this->begin(), this->end(), _obj, [&duals](double sum, const std::pair<int, double>& p) { return sum - p.second * duals[p.first]; });
	}

};




class Constraint : public IloRange {
protected:
	static unsigned int _index;
	unsigned int _id;
	Env* _env;

	NumVar* _SlackVariable;
	double _slackcoef;

	Solver::ConstraintType _type;

	void init() { addVars(_env->getVars()); }

	bool _ended;
public:

	Constraint(Env& env, double lb, double ub, Solver::ConstraintType t) : 
		IloRange(static_cast<IloEnv&>(env), lb, ub), _env(&env), _type(t), _id(_index++),_SlackVariable(nullptr), _slackcoef(0), _ended(false) {
		string name = "C" + to_string(_id);
		this->setName(name.c_str());
	}

	Constraint(Env& env, double lb, double ub, Solver::ConstraintType t, string cstname) : 
		IloRange(static_cast<IloEnv&>(env), lb, ub, cstname.c_str()), _env(&env), _type(t), _id(_index++), _SlackVariable(nullptr), _slackcoef(0), _ended(false) {	}

	virtual ~Constraint() { 
		if (!_ended) this->end(); // End the constraint in the environment
	}

	virtual double getCoefficient(int i, int j) const = 0;
	virtual double getCoefficient(NumVar* X) const = 0;
	void addVar(NumVar* X) { this->setLinearCoef(*X, this->getCoefficient(X)); }
	void addVars(vector<NumVar*>& X);
	Solver::ConstraintType getType() const { return _type; }
	unsigned int getId() const { return _id; }
	void setId(unsigned int id) { _id = id; }

	double retrieveCoefficient(NumVar* x);

	pair<NumVar*, double> toStandarForm();

	void setasended() { _ended = true; } // Set the constraint as ended, so it won't be ended again in the destructor

	bool isValidUnderConstraint(const std::vector<NumVar*>& Vars, const std::vector<double>& Vals);

	void print(const std::vector<NumVar*>& Vars, const std::vector<double>& Vals);

	NumVar* getSlackVariable() { return _SlackVariable; }
	double getSlackCoefficient() const { return _slackcoef; }

};

class ConstraintVarBound : public Constraint {
	NumVar* _X;

public:

	// Only one constructor with an additional parameter to specify if it's an upper or lower bound
	ConstraintVarBound(Env& env, NumVar* X, double bound, bool isUpper)
		: Constraint(env, isUpper ? -IloInfinity : bound, isUpper ? bound : IloInfinity, isUpper ? Solver::OnEdge : Solver::OnEdge), _X(X) {
		string name = (isUpper ? "UB" : "LB") + to_string(X->getId());
		this->setName(name.c_str());
		init();
	}
		
	~ConstraintVarBound() {}
	double getCoefficient(int i, int j) const { return _X->geti() == i && _X->getj() == j ? 1.0 : 0.0; }
	double getCoefficient(NumVar* X) const {
		if (X == _SlackVariable) return _slackcoef;
		return (X == _X ? 1.0 : 0.0); 
	}
};

//class ConstraintNodeDegree : public Constraint {
//	int _i;	//Node i
//	int _degree;
//public:
//	ConstraintNodeDegree(Env& env, int i, int degree=2) : Constraint(env, degree, degree, Solver::OnNode, "D" + to_string(i)), _i(i), _degree(degree) {}
//
//	~ConstraintNodeDegree() {}
//	double getCoefficient(int i, int j) const { return (i == _i || j == _i); }
//	double getCoefficient(NumVar* X) const { 
//		if (X == _SlackVariable) return _slackcoef;
//		return X->getType()== Solver::EdgeVar && (X->geti() == _i || X->getj() == _i); 
//	}
//};

class ConstraintNodeDegreeOut : public Constraint {
	int _i;	//Node i
	int _degree;
public:
	ConstraintNodeDegreeOut(Env& env, int i, int degree = 1) : Constraint(env, degree, degree, Solver::OnNode, "DOut" + to_string(i)), _i(i), _degree(degree) {
		init();
	}

	~ConstraintNodeDegreeOut() {}
	double getCoefficient(int i, int j) const { return (i == _i); }
	double getCoefficient(NumVar* X) const { 
		if (X == _SlackVariable) return _slackcoef;
		return X->getType() == Solver::EdgeVar && X->geti() == _i;
	}
};

class ConstraintNodeDegreeIn : public Constraint {
	int _j;	//Node i
	int _degree;
public:
	ConstraintNodeDegreeIn(Env& env, int j, int degree = 1) : Constraint(env, degree, degree, Solver::OnNode, "DIn" + to_string(j)), _j(j), _degree(degree) {
		init();
	}

	~ConstraintNodeDegreeIn() {}
	double getCoefficient(int i, int j) const { return (j == _j); }
	double getCoefficient(NumVar* X) const { 
		if (X == _SlackVariable) return _slackcoef;
		return X->getType() == Solver::EdgeVar && X->getj() == _j; }
};

class ConstraintLoadPropagation : public Constraint {
	int _i;	//Node i
	int _j;	//Node j
public:
	ConstraintLoadPropagation(Env& env, int i, int j) : Constraint(env, -IloInfinity, env.getInstance()->VehicleCap(), Solver::OnEdge, "L" + to_string(i) + "_" + to_string(j)), _i(i), _j(j) {
		//init();
		//manually add the variables to the constraint: faster
		
		NumVar* xij = _env->getEdgeVar(_i, _j);
		NumVar* ui = _env->getCapVar(_i);
		NumVar* uj = _env->getCapVar(_j);
		if(xij==NULL || ui==NULL || uj==NULL) 
			throw string("Error: Edge or Capacity variable not found for nodes " + to_string(_i) + " and " + to_string(_j));
		addVar(xij);
		addVar(ui);
		addVar(uj);

		//vector<NumVar*>& vars = _env->getVars();
		//for (auto& X : vars) {
		//	if (X->getType() == Solver::EdgeVar && (X->geti() == _i || X->getj() == _j)) {
		//		this->addVar(X);
		//	}
		//	else if (X->getType() == Solver::CapVar && (X->geti() == _i || X->geti() == _j)) {
		//		this->addVar(X);
		//	}
		//}
		//addVars(vars); // Add all variables to the constraint
	}
	~ConstraintLoadPropagation() {}
	double getCoefficient(int i, int j) const {
		if (i == _i && j == _j) return _env->getInstance()->VehicleCap() + _env->getInstance()->Node(j)->demand();
		else return 0.0;
	}
	double getCoefficient(NumVar* X) const { 
		if (X == _SlackVariable) return _slackcoef;
		Solver::VarType type = X->getType();
		if (type == Solver::EdgeVar) return getCoefficient(X->geti(), X->getj());
		else if (type == Solver::CapVar) {
			if (X->geti() == _i) return 1.0;
			else if (X->geti() == _j) return -1.0;
			else return 0.0;
			
		}
		else return 0.0;
	}
};

class ConstraintSubtourElimination : public Constraint {
	unordered_set<int> _S;
public:
	ConstraintSubtourElimination(Env& env, const vector<int>& S) : Constraint(env, -IloInfinity, int(S.size() - 1), Solver::OnEdge, "SC_" + to_string(_id)),
		_S(S.begin(), S.end()) {

		init();

		Instance* I = _env->getInstance();
		double vcap = I->VehicleCap();
        double demand = accumulate(_S.begin(), _S.end(), 0.0, [&](double sum, int i) {  
            return sum + I->Node(i)->demand();  
        });
		double r = ceil(demand / vcap);

		this->setUB(S.size() - 1);
	}

	~ConstraintSubtourElimination() {}
	double getCoefficient(int i, int j) const {
		if (_S.find(i) != _S.end() && _S.find(j) != _S.end()) return 1.0;
		return 0.0;
	}
	double getCoefficient(NumVar* X) const {
		if (X == _SlackVariable) return _slackcoef;
		else if (X->getType() == Solver::EdgeVar) return getCoefficient(X->geti(), X->getj());
		else return 0.0;
	}

};

class ConstraintConnectivity : public Constraint {
public:
	enum class direction { Out, In };
private:	

	unordered_set<int> _S;
	direction _dir;
public:

	ConstraintConnectivity(Env& env, const vector<int>& S, direction dir) : Constraint(env, -IloInfinity, IloInfinity, Solver::OnEdge, "Co" + to_string(_id)),
		_S(S.begin(), S.end()), _dir(dir) {


		Instance* I = _env->getInstance();
		double vcap = I->VehicleCap();
		double demand = 0.0;
		for(int i : _S) {
			if (I->Node(i)->demand() < 0) throw string("Error: Negative demand for node " + to_string(i));
			demand += I->Node(i)->demand();
		}
		//assert(vcap > 0 && "Vehicle capacity must be positive");
		//if (demand <= 0) throw string("Error: Negative demand for the set of nodes in the constraint");
		double r = ceil(demand / vcap);

		this->setLB(r);

		init();

	}

	~ConstraintConnectivity() {}
	double getCoefficient(int i, int j) const {
		if (_dir==direction::Out && _S.find(i) != _S.end() && _S.find(j) == _S.end()) return 1.0;
		else if(_dir == direction::In && _S.find(i) == _S.end() && _S.find(j) != _S.end()) return 1.0;
		return 0.0;
	}
	double getCoefficient(NumVar* X) const {
		if (X == _SlackVariable) return _slackcoef;
		else if (X->getType() == Solver::EdgeVar) return getCoefficient(X->geti(), X->getj());
		else return 0.0;
	}

};


//class ConstraintInfPath : public Constraint {
//	map<int, int> _S;
//	string name() const { return "P_" + to_string(_id); }
//public:
//	ConstraintInfPath(Env& env, const vector<int>& S) : Constraint(env, 0, S.size() - 1, Env::OnEdge) {
//		for (int i = 0; i < S.size(); i++) {
//			_S[S[i]] = S[(i + 1) % S.size()];
//		}
//	}
//
//	~ConstraintInfPath() {}
//	double getCoefficient(int i, int j) const {
//		auto iti = _S.find(i);
//		auto itj = _S.find(j);
//		if (iti != _S.end() && itj != _S.end() && (iti->second == j || itj->second == i)) return 1.0;
//		return 0.0;
//	}
//};   
//
////class ConstraintMatching : public Constraint {
////	unordered_set<int> _S;
////	vector<Edge> _T;
////	string name() const { return "M_" + to_string(_id); }
////public:
////	ConstraintMatching(Env& env, const vector<int>& S, const vector<Edge>& T) : Constraint(env, -IloInfinity, S.size() + (double)(T.size() - 1) / 2.0, Env::OnEdge), _S(S.begin(), S.end()), _T(T) {
////		for (auto& e : _T) if (e.i < e.j) swap(e.i, e.j);
////	}
////	~ConstraintMatching() {}
////	double getCoefficient(int i, int j) const {
////		//return (_S.find(i) != _S.end() && _S.find(j) != _S.end());
////		if (_S.find(i) != _S.end() && _S.find(j) != _S.end()) return 1.0;
////		else {
////			if (i < j) std::swap(i, j);
////			for (auto& e : _T) if (e.i == i && e.j == j) return 1.0;
////			return 0.0;
////		}
////
////	}
////};

class ConstraintGomory : public Constraint {
private:
	string name() const { return "G_" + to_string(_id); }
	vector<Constraint*> _R;
	vector<double> _coeffs;

	mutable std::unordered_map<unsigned int, double> _mapcoef;
public:
	ConstraintGomory(Env& env, double lb, double ub, vector<Constraint*>& R, vector<double>& coeffs, bool autoinit=false) : Constraint(env, lb, ub, Solver::Gomory ), _R(R), _coeffs(coeffs) {
		if (_R.size() != _coeffs.size())
			throw string("Error in the size of the coefficients");

		//We remove from CSTS the constraints with associated coeff 0
		std::vector<Constraint*> new_R;
		std::vector<double> new_coeffs;
		for (size_t i = 0; i < _R.size(); ++i) {
			if (fabs(_coeffs[i]) > EPSILON) {
				new_R.push_back(_R[i]);
				new_coeffs.push_back(_coeffs[i]);
			}
		}
		_R.swap(new_R);
		_coeffs.swap(new_coeffs);


		if(autoinit) init(); //Initialize the constraint in the environment

#if DEBUG
		for (auto r : R)
			if (r->getId() >= _id)
				throw string("Error in the order of the Gomory constraints");
#endif
	}


	~ConstraintGomory() {
		_mapcoef.clear(); // Clear the map to free memory
	}

  //  void updateCoef(NumVar* X, double coef) {
		////_mapcoef[X] = coef; // Store the coefficient in the map for fast access later
  //      // Only update if key exists, otherwise create key
		////unsigned int id = X->getId();
  ////      auto it = _mapcoef.find(id);
  ////      if (it != _mapcoef.end()) {
  ////          it->second = coef; // Update existing key
  ////      } else {
  ////          _mapcoef.emplace(id, coef); // Create new key
  ////      }
  //      this->setLinearCoef(*X, coef); // Call the base class method to set the coefficient in the model
  //      
  //  }

	double getCoefficient(int i, int j) const {
		double sum = 0.0;
		size_t n = _R.size();
		for (size_t k = 0; k < n; k++) {
			sum += _coeffs[k] * _R[k]->getCoefficient(i, j);

		}
		//return sum;
		return sum;//ORUtils::FractionalPart(sum);
	}

	double getCoefficient(NumVar* X) const {
		//throw string("Error: getCoefficient(NumVar*) not implemented for ConstraintGomory. Use retrieveCoefficient instead.");

		//Improve speed at the cost of memory used
		unsigned int id = X->getId();
		auto it = _mapcoef.find(id);
		if (it != _mapcoef.end()) return it->second;


		if (X == _SlackVariable) return _slackcoef;
		double sum = 0.0;
		size_t n = _R.size();
		for (size_t k = 0; k < n; k++) {
			sum += _coeffs[k] * _R[k]->getCoefficient(X);

		}
		sum = ORUtils::FractionalPart(sum);

		_mapcoef[id] = sum;

		return sum;


	}


};

