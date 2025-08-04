#pragma once


#include "../LPSolvers/LPSolver.hpp"
#include <queue>
#include <vector>

typedef std::map<NumVar*, Solver::VarStatus> VarStatusMap;

// Class B&B Node
class BBNode{
public:
	

private:
	static unsigned int _index;
	unsigned int _id;			//Node ID
	unsigned int _depth;			//Depth of the Node
	std::map<NumVar*, Solver::VarStatus> _statusMap;

	Solver::NodeStatus _status;				//Node Status

	double _objval;					//Objective Value
	BBNode* _parent;				//Parent Node
	BBNode* _left;					//Down Node
	BBNode* _right;					//Up Node



public:

	static void applyVarStatus(NumVar* X, Solver::VarStatus status);
	static void revokeVarStatus(NumVar* X);

	static void applyNodeVarStatus(BBNode* node);
	static void revokeNodeVarStatus(BBNode* node);


	//Constructor for the Root Node
	BBNode(double objvalue = DBL_MAX) : _parent(NULL), _status(Solver::NodeStatus::PendingNode),
		_objval(objvalue), _left(NULL), _right(NULL), _depth(0), _id(_index++) {}

	//Constructor for a Node
	BBNode(NumVar* X, Solver::VarStatus varstatus, BBNode* parent, double objvalue = DBL_MAX) :
		_parent(parent), _status(Solver::NodeStatus::PendingNode), _objval(objvalue), _left(NULL), _right(NULL), _id(_index++) {
		_depth = (parent == NULL) ? 0 : parent->_depth + 1;
		_statusMap[X] = varstatus;
	}

	//Constructor for a Node with VarStatusMap
	BBNode(VarStatusMap& statusMap, BBNode* parent, double objvalue = DBL_MAX) :
		_parent(parent), _status(Solver::NodeStatus::PendingNode), _objval(objvalue), _left(NULL), _right(NULL), _id(_index++) {
		_depth = (parent == NULL) ? 0 : parent->_depth + 1;
		_statusMap = statusMap;
	}

	~BBNode();
	unsigned int getId() const { return _id; }
	unsigned int getDepth() const { return _depth; }

	Solver::NodeStatus getStatus() const { return _status; }
	Solver::VarStatus getNodeDir() const { return _statusMap.begin()->second; }
	double getObjVal() const { return _objval; }
	BBNode* getParent() const { return _parent; }
	BBNode* getLeft() const { return _left; }
	BBNode* getRight() const { return _right; }
	void setStatus(Solver::NodeStatus status) { _status = status; }
	void setObjVal(double objval) { _objval = objval; }
	void generateChildren(NumVar* X);
	void generateChildren(VarStatusMap& statusMapLeft, VarStatusMap& statusMapRight);

	void setStatus(NumVar* X, Solver::VarStatus status);

	//Ordering operator
	bool operator<(const BBNode& node) const { return _objval < node._objval; }

	void setLocalVarBound();

	bool isRoot() const { return _parent == NULL; }
	bool isLeaf() const { return _left == NULL && _right == NULL; }

};

struct BBNodePtrCompare {
    bool operator()(const BBNode* lhs, const BBNode* rhs) const {
        // Min-heap: node with smallest objval (best) is on top
        return lhs->getObjVal() > rhs->getObjVal();
    }
};

class IPSolver : public Solver {

	LPSolver* _LP;
	void solvemethod(Solution* S) {
		_LP->convertToIP();
		_LP->setparam(Solver::TimLim, _timlim);
		_LP->setparam(Solver::Threads, _threads);
		_LP->setparam(Solver::Gap, _mingap);
		_LP->enableOutput(true);
		_LP->solve(S);
	}



public:
	IPSolver(LPSolver* LP) : Solver(LP->instance(), "IP"), _LP(LP) { _LP->convertToIP(); }
	~IPSolver() {}
	
	double gap() { return _LP->gap(); }
	Solution* recoversolution() { return _LP->recoversolution(); }


};

class BBSolver : public Solver {

	LPSolver* _LP;
	std::priority_queue<BBNode*, std::vector<BBNode*>, BBNodePtrCompare> _BBNodes;
	Solution* _bestsol;
	double _bestUB;

	void solvemethod(Solution* S);

	bool updateBestTour(Solution* S);

	void addNode(BBNode* node);

	vector<FlowEdge*> getBranchingCandidates(Solution* S);

	void printlog(int freq, BBNode* N=NULL);

	std::pair< std::vector<NumVar*>, std::vector<NumVar*> > _splitVars(Solution* sol, int i) const;


public:
	BBSolver(LPSolver* LP) : Solver(LP->instance(), "BB"), _LP(LP), _bestsol(NULL), _bestUB(DBL_MAX) { }
	~BBSolver() {
        while (!_BBNodes.empty()) { delete _BBNodes.top(); _BBNodes.pop(); }
    }
	double gap(){
		if (_bestsol == NULL) return 1.0;
		else if (_BBNodes.empty()) return 0.0;
		else return ( _bestsol->fitness() - _BBNodes.top()->getObjVal() ) / _bestsol->fitness();
	}
	Solution* recoversolution();

    // Greedy repair using NodeSol's flows
    Solution* greedyConstructSolution(Solution* NodeSol);
};                                               