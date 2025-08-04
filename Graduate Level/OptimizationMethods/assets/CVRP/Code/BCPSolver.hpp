#pragma once



#include "Solver.hpp"
#include "ColumnGen.hpp"

class BCPSolver;

//Class B&B Node
class BBNode
{
private:
	static unsigned int _index;
	unsigned int _id;			//Node ID
	unsigned int _depth;			//Depth of the Node
	CGNumVar* _X;					//Branching Variable
	CGEnv::NodeStatus _status;
	CGEnv::VarStatus _varstatus;
	double _objval;					//Objective Value
	BBNode* _parent;				//Parent Node
	BBNode* _left;					//Down Node
	BBNode* _right;					//Up Node

	vector<CGNumVar*> _prunedVars;

public:

	//Constructor for the Root Node
	BBNode(CGEnv& env, double objvalue = DBL_MAX) : _X(NULL), _parent(NULL), _status(CGEnv::PendingNode), _varstatus(CGEnv::Free),
		_objval(objvalue), _left(NULL), _right(NULL), _depth(0), _id(_index++) {
		_X = NULL;
	}

	//Constructor for a Node
	BBNode(CGNumVar* X, CGEnv::VarStatus varstatus, BBNode* parent, double objvalue = DBL_MAX) : _X(X),
		_parent(parent), _status(CGEnv::PendingNode), _varstatus(varstatus), _objval(objvalue), _left(NULL), _right(NULL), _id(_index++) {
		_depth = (parent == NULL) ? 0 : parent->_depth + 1;
	}
	~BBNode() {
		//Delete Children
		if (_left != NULL)
			delete _left;
		if (_right != NULL)
			delete _right;

		//Delete Parent
		if (_parent != NULL) {
			if (_parent->_left == this) _parent->_left = NULL;
			if (_parent->_right == this) _parent->_right = NULL;
			if (_parent->_left == NULL && _parent->_right == NULL) delete _parent;
		}

	}
	unsigned int getId() const { return _id; }
	unsigned int getDepth() const { return _depth; }
	CGNumVar* getVar() const { return _X; }
	CGEnv::NodeStatus getStatus() const { return _status; }
	CGEnv::VarStatus getNodeDir() const { return _varstatus; }
	double getObjVal() const { return _objval; }
	BBNode* getParent() const { return _parent; }
	BBNode* getLeft() const { return _left; }
	BBNode* getRight() const { return _right; }
	void setStatus(CGEnv::NodeStatus status) { _status = status; }
	void setObjVal(double objval) { _objval = objval; }
	void generateChildren(CGNumVar* X) {
		_left = new BBNode(X, CGEnv::Down, this, _objval);
		_right = new BBNode(X, CGEnv::Up, this, _objval);
	}

	void setAsPruned(CGNumVar* x) {
		_prunedVars.push_back(x);
	}

	//Ordering operator
	bool operator<(const BBNode& node) const {
		return _objval < node._objval;
	}

	void setVarBound() {
		if (_varstatus == CGEnv::Down) {
			_X->setLB(0);
			_X->setUB(0);
		}
		else if (_varstatus == CGEnv::Up) {
			_X->setLB(1);
			_X->setUB(1);
		}

		for (auto x : _prunedVars) {
			x->setLB(0);
			x->setUB(0);
		}
	}

	bool isRoot() const { return _parent == NULL; }
	bool isLeaf() const { return _left == NULL && _right == NULL; }

};


class CGMasterProblem
{
	CGEnv _env;
	IloModel _model;
	IloCplex _cplex;
	IloObjective _obj;
	vector< CGNumVar* > _X;
	vector< CGConstraint* > _constraints;
	int _maxColPerRound;

	//We store the columns in a map
	unsigned int _ncols;
	CGNumColumn** _columns;

	bool isEdgeVarInModel(int i, int j) {
		return getCol(i, j)->isInModel();
	}

	double getObjCoeff(int i, int j) { 
		Instance* I = _env.getInstance();
		return I->travellingtime(I->Node(i), I->Node(j)); 
	}

	void setBounds(BBNode* N);

	CGNumColumn* getCol(int i, int j) {
		if (i == j) return NULL;
		CGNumColumn* col = _columns[getColIndex(i, j)];
		assert(col != NULL);
		assert((col->geti() == i && col->getj() == j) || (col->geti() == j && col->getj() == i));
		return _columns[getColIndex(i, j)];

	}

	int getColIndex(int i, int j) {
		if (i < j) swap(i, j);
		return (i * (i - 1)) / 2 + j;
	}

	vector<CGNumVar*> pruneNode(BBNode* N);

public:
	CGMasterProblem(Instance* I);
	~CGMasterProblem();
	void addVar(int i, int j);
	void addConstraint(CGConstraint* constraint);
	void addConstraints(vector<CGConstraint*>& constraints);
	void setupConstraint(CGConstraint* constraint) { constraint->addVars(_X); }
	CGEnv::SolverStatus solve(BBNode* N, double UB=DBL_MAX);
	double getObjVal() const { return _cplex.getObjValue(); }
	double getVarVal(CGNumVar& X) const { return _cplex.getValue(X); }
	IloNumArray getVarVals();
	IloNumArray getDuals();
	void setMaxColPerRound(int max) { _maxColPerRound = max; }
	CGEnv& getEnv(){ return _env; }
	int getNumVars() const { return int(_X.size()); }
	int getNumConstraints() const { return int(_constraints.size()); }
	CGConstraint* getConstraint(int i) { return _constraints[i]; }
	CGNumVar* getVar(int i) { return _X[i]; }

	bool isSlackVar(int k) { return _X[k]->getType() == CGEnv::SlackVar; }
	bool isEdgeVar(int k) { return _X[k]->getType() == CGEnv::EdgeVar; }

	vector<double> getInvRow(int k);
	vector<double> getInvBasisRow(int k);
	double getInvRHS(int k);
	pair<int*,double*> getBHead();
	IloCplex::BasisStatusArray getStatus();

	//We compute the reduced cost of a variable
	double getReducedCost(int i, int j, IloNumArray& duals);
	double getReducedCost(CGNumColumn* col, IloNumArray& duals);
	void addCol(CGNumColumn* col);

	//We compute the distance to a cut
	double getDistanceToCut(CGConstraint& cut, IloNumArray& vals);

	//Criteria for Cut acceptance
	bool acceptCut(CGConstraint& cut, IloNumArray& vals) {
		return getDistanceToCut(cut, vals) > CGEnv::DELTA;
	}

	double getLength(int i, int j) { return getCol(i, j)->getObj(); }
};


class BBFlowGraph;

class FlowEdge : public unordered_set<int> {
	friend class BBFlowGraph;
	int id;			//Edge ID
	int i;			//Node i
	int j;			//Node j	
	FlowEdge* r;	//Reverse Edge
	double w;		//Weight
	double f;		//Flow
public:
	FlowEdge() : id(-1), i(-1), j(-1), w(0), f(0), r(NULL) {}
	~FlowEdge() {}

	bool contains(int k) { return find(k) != end(); }
};

//Supported Grapg from LP solution for max flow
class BBFlowGraph {

	vector<FlowEdge*> _edges;
	vector<vector<int>> _adj;	
	vector<int> _pred;
	unordered_set<int> _sources;
	unordered_set<int> _sinks;
	int _nnodes;
	int _nedges;

	bool tryMergeEdge(int i, int j, double wij, double wji);

public:

	BBFlowGraph(int nnodes) : _nnodes(nnodes), _nedges(0) {
		_adj.resize(nnodes);
		_pred.resize(nnodes, -1);
	}

	~BBFlowGraph() {
		for (auto e : _edges) delete e;
	}

	double residualflow(FlowEdge* e) { return e->w + e->f; }

	double getFlow(const vector<int>& shore1, const vector<int>& shore2);

	double getLeavingFlow(const vector<int>& shore);

	void addEdge(int i, int j, double wij, double wji) {
		FlowEdge* eij = new FlowEdge();
		FlowEdge* eji = new FlowEdge();
		eij->id = _nedges++;
		eij->i = i;
		eij->j = j;
		eij->w = wij;
		eij->f = 0;
		eji->id = _nedges++;
		eji->i = j;
		eji->j = i;
		eji->w = wji;
		eji->f = 0;
		eij->r = eji;
		eji->r = eij;
		_edges.push_back(eij);
		_edges.push_back(eji);
		_adj[i].push_back(eij->id);
		_adj[j].push_back(eji->id);
	}

	void addSource(int i) { _sources.insert(i); }
	void addSink(int i) { _sinks.insert(i); }
	void resetSources() { _sources.clear(); }
	void resetSinks() { _sinks.clear(); }

	void resetFlow() { 
		for (auto e : _edges) e->f = 0; 
		fill(_pred.begin(), _pred.end(), -1);
	}

	double getFlow(vector<FlowEdge*>& P);

	vector<FlowEdge*> getPath();

	double maxFlow();

	double getShores(vector<int>& shore1, vector<int>& shore2, bool strict = false);

	void exportGraph(string name);
};

//Class Cutting Plane Solver
class BBCuttingPlaneSolver {
	CGMasterProblem* _master;
	Instance* _I;

    // Add a pointer to function for heuristic
	std::function<bool(IloNumArray&)> heuristicFunction;

	int findparent(int i, vector<int>& parent) {
		if (parent[i] == -1) return i;
		else return parent[i] = findparent(parent[i], parent);
	}
	void unite(int i, int j, vector<int>& parent, vector<int>& rank) {
		int s1 = findparent(i, parent);
		int s2 = findparent(j, parent);
		if (s1 != s2) {
			if (rank[s1] < rank[s2])  parent[s1] = s2;
			else if (rank[s1] > rank[s2]) parent[s2] = s1;
			else {
				parent[s2] = s1;
				rank[s1] += 1;
			}
		}
	}

public:

	BBCuttingPlaneSolver(CGMasterProblem* M, std::function<bool(IloNumArray&)> heuristic = 0) : _I(M->getEnv().getInstance()),
		_master(M), heuristicFunction(heuristic) {}
	~BBCuttingPlaneSolver() {};

	bool getFlowCuts(IloNumArray& vals, vector<CGConstraint*>& Cuts);

	bool getSECCuts(IloNumArray& vals, vector<CGConstraint*>& Cuts);

	bool getGomoryCuts(IloNumArray& vals, vector<CGConstraint*>& Cuts);

	bool getMatchingCut(IloNumArray& vals, vector<CGConstraint*>& Cuts);

	bool getChavtalCuts(IloNumArray& vals, vector<CGConstraint*>& Cuts);

	bool getInfeasiblePathCuts(IloNumArray& vals, vector<CGConstraint*>& Cuts, double IncObj);

	void solve(BBNode* N);
};





//Class B&B Solver
class BCPSolver : public Solver
{
	CGMasterProblem* _master;
	BBCuttingPlaneSolver* _cuttingplanesolver;
	vector<BBNode*> _BBNodes;
	vector<int> _besttour;
	double _bestobjval;
	
	void solvemethod(Solution* S);

	void initMasterProblem(const int maxGenTours=10, const int k = 3);

	void initMasterProblemRandom(const int maxGenTours = 2);

	void addAllVars();	//For debugging purposes

	bool updateBestTour(IloNumArray& vals);

	void addNode(BBNode* node);

	bool applyHeuristic(IloNumArray& vals);

	vector<int> extractTour(IloNumArray& vals);

	bool addFeasibleTour(vector<int>& tour);

	vector<CGNumVar*> getBranchingCandidates(IloNumArray& vals);	

public:

	BCPSolver(Instance* I);
	~BCPSolver();
	double gap() { 
		if (_besttour.empty()) return 1.0;
		else if (_BBNodes.empty()) return 0.0;
		else {
			double bestbound = _BBNodes.back()->getObjVal();
			return (_bestobjval - bestbound) / bestbound;
		}
	}
	Solution* recoversolution();



};

