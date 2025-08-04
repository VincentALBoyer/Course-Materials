#pragma once

//This file is the main header file for the BBColumnGen library
//It includes all the classes and functions that are part of the column generation library
//It is build over the CPLEX library

#include "Instance.hpp"

struct CGEdge {
	int i;
	int j;
	double w;

	bool operator==(const CGEdge& other) const {
		return i == other.i && j == other.j;
	}

	// Define the order relationship
	bool operator<(const CGEdge& other) const {
		if (i != other.i) return i < other.i;
		if (j != other.j) return j < other.j;
		return w < other.w;
	}
};

class CGEnv : public IloEnv {
	Instance* _I;

	static double FLOOR(double x) {
		double val = floor(x);
		if (fabs(x - val) > 1.0 - EPSILON) {
			return val + 1;
		}
		else return val;
	}

public:
	enum VarStatus { Free, Up, Down };
	enum SolverStatus { Optimal, Feasible, Infeasible, Unknown };
	enum ConstraintType { OnNode, OnEdge, Gomory, NotSpecified };
	enum VarType { EdgeVar, SlackVar };
	enum NodeStatus { InfeasibleNode, FeasibleNode, IntegerNode, IntegerInfeasibleNode, PendingNode, LeafNode };


	//Epsilon for double comparison
	const static double EPSILON;

	//DELTA for Cut Acceptance
	const static double DELTA;

	CGEnv(Instance* I) : IloEnv(), _I(I) {}
	CGEnv() : IloEnv(), _I(NULL) {}
	~CGEnv() { this->end(); }

	void setInstance(Instance* I) { _I = I; }
	Instance* getInstance() const { return _I; }
	static bool isInteger(double x) { return fabs(x - round(x)) < EPSILON; }
	static bool isFractional(double x) { return !isInteger(x); }
	static bool isZero(double x) { return fabs(x) < EPSILON; }
	//static bool fractionalPart(double x) { return x - floor(x); }
	static bool isInteger(IloNumArray& x) {
		for (int i = 0; i < x.getSize(); i++) {
			if (!isInteger(x[i])) return false;
		}
		return true;
	}

	static double FractionalPart(double x) { return x - FLOOR(x); }
	static double PositiveFractionalPart(double x) {
		//if (x < 0) return  x - ceil(x);
		return x - FLOOR(x);
		
	}

	
};




class CGNumVar : public IloNumVar {

protected:
	static unsigned int _index;

	unsigned int _id;

	int _i;
	int _j;

	void consolidate() {
		if (_type == CGEnv::EdgeVar){
			if (_i < _j) swap(_i, _j);
			setName(CGNumVar::nameEdgeGen(_i, _j).c_str());
		}
		else {
			setName(CGNumVar::nameSlackGen().c_str());
			_i = _j = -1;	
		}
	}

	CGEnv::VarType _type;

	CGEnv* _env;

public:

	static string nameSlackGen() {
		return "S_" + to_string(_index);
	}

	static string nameEdgeGen(int i, int j) {
		int ii = max(i, j);
		int jj = min(i, j);
		return "X" + to_string(ii) + "_" + to_string(jj);
	}

	//static unsigned int hash(int i, int j) { return i * 1000 + j; }

	CGNumVar(CGEnv& env) : IloNumVar(static_cast<IloEnv&>(env), 0, IloInfinity), _id(_index++), _env(&env), _type(CGEnv::SlackVar) {
		consolidate();
	}
	CGNumVar(CGEnv& env, int i, int j) : IloNumVar(static_cast<IloEnv&>(env), 0, 1), _id(_index++), _env(&env), _type(CGEnv::EdgeVar), _i(i), _j(j) {
		consolidate();
	}
	CGNumVar(CGNumVar& var) : IloNumVar(var), _id(var._id), _env(var._env), _type(var._type), _i(var._i), _j(var._j) {}
	//CGNumVar(CGEnv& env, IloNumColumn& col) : IloNumVar(col, 0, IloInfinity), _id(_index++), _env(&env), _type(CGEnv::SlackVar) {
	//	consolidate();
	//}
	CGNumVar(CGEnv& env, IloNumColumn& col, int i, int j) : IloNumVar(col, 0, 1), _id(_index++), _env(&env), _type(CGEnv::EdgeVar), _i(i), _j(j) {
		consolidate();
	}
	~CGNumVar() { this->end(); }

	unsigned int getId() const { return _id; }

	CGEnv::VarType getType() const { return _type; }	

	void apply(CGEnv::VarStatus status) {
		if (status == CGEnv::Up) {
			this->setLB(1);
			this->setUB(1);
		}
		else if (status == CGEnv::Down) {
			this->setLB(0);
			this->setUB(0);
		}
		else {
			this->setLB(0);
			this->setUB(1);
		}
	}

	int geti() const { return _i; }
	int getj() const { return _j; }


};

class CGNumColumn : public std::map<int, double>  {

	CGNumVar* _var;

	double _obj;
	int _i;
	int _j;

public:

	CGNumColumn(int i, int j, double obj) : std::map<int, double>(), _i(i), _j(j), _obj(obj), _var(NULL) {
		if (_i < _j) std::swap(_i, _j);
	}
	~CGNumColumn() {}

	double getObj() const { return _obj; }
	int geti() const { return _i; }
	int getj() const { return _j; }
	CGNumVar* getVar() const { return _var; }

	void setObj(double obj) { _obj = obj; }
	void add(int i, double val) {
		if (fabs(val) < CGEnv::EPSILON) return;
		if (isInModel()) return;

		this->insert_or_assign(i, val);
	}

	void setLB(double lb) { _var->setLB(lb); }
	void setUB(double ub) { _var->setUB(ub); }

	bool isInModel() const { return _var != NULL; }
	void setInModel(CGNumVar* x) {
		//this->clear(); //We do not need the information anymore (memory optimization)
		_var = x;
	}

	double getReducedCost(IloNumArray& duals) {
		return std::accumulate(this->begin(), this->end(), _obj, [&duals](double sum, const std::pair<int, double>& p) { return sum - p.second * duals[p.first]; });
	}

};

//class CGNumColumn : public std::vector< pair<int, double> > {
//
//	CGNumVar* _var;
//
//	double _obj;
//	int _i;
//	int _j;
//
//public:
//
//	CGNumColumn(int i, int j, double obj) : std::vector< pair<int, double> >(), _i(i), _j(j), _obj(obj), _var(NULL) {
//		if (_i < _j) std::swap(_i, _j);
//	}
//	~CGNumColumn() {}
//
//	double getObj() const { return _obj; }
//	int geti() const { return _i; }
//	int getj() const { return _j; }
//	CGNumVar* getVar() const { return _var; }
//
//	void setObj(double obj) { _obj = obj; }
//	void add(int i, double val) { 
//		if (fabs(val) < CGEnv::EPSILON) return;
//		if (isInModel()) return;
//		auto it = std::find_if(this->begin(), this->end(), [i](const pair<int, double>& p) { return p.first == i; });
//		if (it == this->end()) {
//			this->emplace_back(i, val);
//		}
//		else if (fabs(val) > CGEnv::EPSILON) {
//			it->second = val;
//		}
//		else {
//			this->erase(it);
//		}
//	}
//	//bool hasKey(int i) {
//	//	for (auto& p : *this) if (p.first == i) return true;
//	//	return false;
//	//}
//	//double getCoefficient(int i) {
//	//	for (auto& p : *this) if (p.first == i) return p.second;
//	//	return 0.0;
//	//}
//
//	bool isInModel() const { return _var != NULL; }
//	void setInModel(CGNumVar* x) {
//		this->clear(); //We do not need the information anymore (memory optimization)
//		_var = x; 
//	}
//
//	double getReducedCost(IloNumArray& duals) {
//		if (_var != NULL) return 0.0;
//		
//		return std::accumulate(this->begin(), this->end(), _obj, [&duals](double sum, pair<int, double> p) { return sum - p.second * duals[p.first]; });
//
//	}
//
//};


class CGConstraint : public IloRange {
protected:
	static unsigned int _index;
	unsigned int _id;
	CGEnv* _env;
	virtual string name() const { return "CGCst" + to_string(_id); }

	CGEnv::ConstraintType _type;

public:

	CGConstraint(CGEnv& env, double lb, double ub, CGEnv::ConstraintType t) : IloRange(static_cast<IloEnv&>(env), lb, ub), _env(&env), _type(t), _id(_index++) {
		this->setName(name().c_str());
	}

	virtual ~CGConstraint() { this->end();  }

	virtual double getCoefficient(int i, int j) const = 0;
	double getCoefficient(CGNumVar& X) const { return getCoefficient(X.geti(), X.getj()); }
	void addVar(CGNumVar& X) { this->setLinearCoef(X, this->getCoefficient(X)); }
	void addVars(vector<CGNumVar*>& X) {
		//for (auto& x : X) this->addVar(*x); 
		IloNumArray coefs(this->getEnv());
		IloNumVarArray vars(this->getEnv());
		for (auto& x : X) {
			coefs.add(this->getCoefficient(*x));
			vars.add(*x);
		}
		this->setLinearCoefs(vars, coefs);
		vars.end();
		coefs.end();
	}
	CGEnv::ConstraintType getType() const { return _type; }
	unsigned int getId() const { return _id; }
	void setId(unsigned int id) { _id = id; }

	double retrieveCoefficient(CGNumVar* x)
	{
		string name = x->getName();
		IloExpr expr = getExpr();
		// Get the linear expression from the range 
		double coefficient = 0.0;
		for (IloExpr::LinearIterator it = expr.getLinearIterator(); it.ok(); ++it) {
			if (it.getVar().getName() == name) {
				coefficient = it.getCoef();
			}
		}
		expr.end(); // Clean up the expression return coefficient;

		return coefficient;
	}

	virtual bool merge(CGConstraint* other) {
		return false;
	}

	virtual void clean() {};
};

class CGConstraintNodeDegree : public CGConstraint {
	int _i;	//Node i
	string name() const { return "D" + to_string(_i); }
public:
	CGConstraintNodeDegree(CGEnv& env, int i) : CGConstraint(env, 2, 2, CGEnv::OnNode), _i(i) {}

	~CGConstraintNodeDegree() {}
	double getCoefficient(int i, int j) const {
		return (i == _i || j == _i);
		//if (i == _i || j == _i) return 1.0;
		//return 0.0;
	}
};

class CGConstraintSubtourElimination : public CGConstraint {
	unordered_set<int> _S;
	string name() const { return "T_" + to_string(_id); }
public:
	CGConstraintSubtourElimination(CGEnv& env, const vector<int>& S) : CGConstraint(env, -IloInfinity, int(S.size() - 1), CGEnv::OnEdge), _S(S.begin(), S.end()) {}

	~CGConstraintSubtourElimination() { }
	double getCoefficient(int i, int j) const {
		//return (_S.find(i) != _S.end() && _S.find(j) != _S.end());
		if (_S.find(i) != _S.end() && _S.find(j) != _S.end()) return 1.0;
		return 0.0;
	}

	bool merge(CGConstraint* other) {
		CGConstraintSubtourElimination* o = dynamic_cast<CGConstraintSubtourElimination*>(other);
		if (o == NULL) return false;
		//We get the  different elements
		vector<int> diff;
		for (auto& s : _S) if (o->_S.find(s) == o->_S.end()) diff.push_back(s);
		for (auto& s : o->_S) if (_S.find(s) == _S.end()) diff.push_back(s);
		if (diff.size() <= 1) {
			_S = o->_S;
			//_S.insert(diff.begin(), diff.end());
			this->setBounds(int(_S.size() - 1), int(_S.size() - 1));
			
			return true;
		}
		else return false;
	}

	void clean() {
		_S.clear();
	}
};

class CGConstraintInfPath : public CGConstraint {
	map<int, int> _S;
	string name() const { return "P_" + to_string(_id); }
public:
	CGConstraintInfPath(CGEnv& env, const vector<int>& S) : CGConstraint(env, 0, S.size()-1, CGEnv::OnEdge) {
		for (int i = 0; i < S.size(); i++) {
			_S[S[i]] = S[(i + 1) % S.size()];
		}
	}

	~CGConstraintInfPath() {  }
	double getCoefficient(int i, int j) const {
		auto iti = _S.find(i);
		auto itj = _S.find(j);
		if (iti != _S.end() && itj != _S.end() && (iti->second == j || itj->second == i)) return 1.0;
		return 0.0;
	}
};

class CGConstraintMatching : public CGConstraint {
	unordered_set<int> _S;
	vector<CGEdge> _T;
	string name() const { return "M_" + to_string(_id); }
public:
	CGConstraintMatching(CGEnv& env, const vector<int>& S, const vector<CGEdge>& T) : CGConstraint(env, -IloInfinity, S.size() + (double)(T.size() - 1) / 2.0, CGEnv::OnEdge), _S(S.begin(), S.end()), _T(T) {
		for (auto& e : _T) if (e.i < e.j) swap(e.i, e.j);
	}
	~CGConstraintMatching() {  }
	double getCoefficient(int i, int j) const {
		//return (_S.find(i) != _S.end() && _S.find(j) != _S.end());
		if (_S.find(i) != _S.end() && _S.find(j) != _S.end()) return 1.0;
		else {
			if (i < j) swap(i, j);
			for (auto& e : _T) if (e.i == i && e.j == j) return 1.0;
			return 0.0;
		}
		
	}
};

class CGConstraintMaxNEdges : public CGConstraint {
	string name() const { return "E_" + to_string(_id); }
public:
	CGConstraintMaxNEdges(CGEnv& env, int n) : CGConstraint(env, n, n, CGEnv::OnEdge) {}

	~CGConstraintMaxNEdges() {  }
	double getCoefficient(int i, int j) const {
		return 1.0;
	}
};

class CGConstraintConnectivity : public CGConstraint {
	unordered_set<int> _S;
	string name() const { return "T_" + to_string(_id); }
public:
	//CGConstraintConnectivity(CGEnv& env, unordered_set<int>& S) : CGConstraint(env, 2, IloInfinity, CGEnv::OnEdge), _S(S) {}
	CGConstraintConnectivity(CGEnv& env, const vector<int>& S) : CGConstraint(env, 2, IloInfinity, CGEnv::OnEdge), _S(S.begin(), S.end()) {}

	~CGConstraintConnectivity() {  }
	double getCoefficient(int i, int j) const {
		bool iInS = _S.find(i) != _S.end();
		bool jInS = _S.find(j) != _S.end();
		if (!iInS && !jInS) return 0.0;
		else if (!iInS || !jInS) return 1.0;
		//if ((!iInS && jInS)|| (iInS && !jInS)) return 1.0;
		return 0.0;
	}
};

class CGConstraintGomory : public CGConstraint {
private:
	string name() const { return "G_" + to_string(_id); }
	vector<CGConstraint*> _R;
	vector<double> _coeffs;
public:
	CGConstraintGomory(CGEnv& env, double lb, double ub, vector<CGConstraint*>& R, vector<double>& coeffs) : CGConstraint(env, lb, ub, CGEnv::Gomory), _R(R), _coeffs(coeffs) {
		if (_R.size() != _coeffs.size())
			throw string("Error in the size of the coefficients");

		//We remove from CSTS the constraints with associated coeff 0
		auto it = _R.begin();
		auto itc = _coeffs.begin();
		while (it != _R.end()) {
			if (fabs(*itc) > CGEnv::EPSILON) {
				it++;
				itc++;
			}
			else {
				it = _R.erase(it);
				itc = _coeffs.erase(itc);
			}
		}

		_R.shrink_to_fit();
		_coeffs.shrink_to_fit();

#if DEBUG
		for (auto r : R)
			if (r->getId() >= _id)
				throw string("Error in the order of the Gomory constraints");
#endif
	}
	
	
	~CGConstraintGomory() {  }
	double getCoefficient(int i, int j) const {
		double sum = 0.0;
		size_t n = _R.size();
		for (size_t k = 0; k < n; k++) {
			sum += _coeffs[k] * _R[k]->getCoefficient(i, j);
			
		}
		//return sum;
		return CGEnv::FractionalPart(sum);
	}

	double getCoefficient(CGNumColumn* C) const {
		double sum = 0.0;
		size_t n = _R.size();
		//for (size_t k = 0; k < n; k++) {
		//	sum += _coeffs[k] * C->getCoefficient( _R[k]->getId());
		//}
		
		auto it = C->begin();
		for (size_t k = 0; k < n; k++) {
			int cstid = _R[k]->getId();
			while (it != C->end() && it->first < cstid) it++;
			if (it != C->end() && it->first == cstid) {
				sum += _coeffs[k] * it->second;
			}
			if (it == C->end()) break;
		}

		return CGEnv::FractionalPart(sum);
	}

	void clean() {
		_R.clear();
		_coeffs.clear();
	}
	
	
};

class CGConstraintChvatal : public CGConstraint {
private:
	string name() const { return "C_" + to_string(_id); }
	vector<CGConstraint*> _R;
	vector<double> _coeffs;
public:
	//CGConstraintGomory(CGEnv& env, double lb, double ub) : CGConstraint(env, lb, ub, CGEnv::NotSpecified), _id(_index++) {}
	CGConstraintChvatal(CGEnv& env, double ub, vector<CGConstraint*>& R, vector<double>& coeffs) : CGConstraint(env, -IloInfinity, floor(ub), CGEnv::OnEdge), _R(R), _coeffs(coeffs) {}
	~CGConstraintChvatal() {  }
	double getCoefficient(int i, int j) const{
		//We do the dot product of the coefficients
		double sum = 0.0;
		size_t n = _R.size();
		for (size_t k = 0; k < n; k++) {
			sum += _coeffs[k] * _R[k]->getCoefficient(i, j);
		}
		//return sum;
		return floor(sum);
	}
};

