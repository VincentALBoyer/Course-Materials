#pragma once

#include "Instance.hpp"


class Solution;
class FlowEdge  {
    friend class Solution;
    int id;			//Edge ID
    int i;			//Node i
    int j;			//Node j	
    FlowEdge* r;	//Reverse Edge
	double rf;		//edge residual flow
    double f;		//Flow
public:
    FlowEdge() : id(-1), i(-1), j(-1), rf(0), f(0), r(NULL) {}
	FlowEdge(int id, int i, int j, double f) : id(id), i(i), j(j), rf(0), f(f), r(NULL) {}
    ~FlowEdge() {}

	int getId() const { return id; }
	int geti() const { return i; }
	int getj() const { return j; }  
	double getFlow() const { return f; }

};


class Solution : protected vector<FlowEdge*> {
    
    Instance* _I;

	vector<vector<int>> _adj;  //Note: maybe change to vector<unordered_set<int>> for better performance
    vector<int> _pred;

    vector<FlowEdge*> getPath(const vector<int>& sources, const vector<int>& sinks);
    double residualflow(FlowEdge* e) { return e->rf; }
    void resetResidualFlow();
	std::vector<int> reachableNodes(const std::vector<int>& sources);

public:
    Solution(Instance* I);;
    Solution(string filename, Instance* I);
	Solution(Solution* S);
    ~Solution();

    int addEdge(int i, int j, double fij , double fji);
	int addEdge(int i, int j, double f = 1.0) { return addEdge(i, j, f, 0.0); }

	double getFlow(int i, int j);

    double getFlow(const std::vector<int>& sources, const std::vector<int>& sinks);

	vector<FlowEdge*> getAdjacentEdges(int i) const;
    
    double fitness();

	bool isInteger();

    vector< vector<int> > extractroutes();
    
    bool isFeasible(bool displayerror = true);
    
    void print();

	void exportGraph(const std::string& filename="Graph.txt");

	static Solution* readSolution(const std::string& filename, Instance* I);

    
    // Compute max flow from sources to sinks, returns pair of shores (reachable, not reachable)
    std::pair<std::vector<int>, std::vector<int>> maxFlow(const std::vector<int>& sources, const std::vector<int>& sinks);

    // Expose iterator and const_iterator from std::vector<int>
    using std::vector<FlowEdge*>::begin;
    using std::vector<FlowEdge*>::end;
    using std::vector<FlowEdge*>::cbegin;
    using std::vector<FlowEdge*>::cend;

    // Expose size and other const functions
    using std::vector<FlowEdge*>::size;
    using std::vector<FlowEdge*>::operator[];
    using std::vector<FlowEdge*>::at;
    using std::vector<FlowEdge*>::front;
    using std::vector<FlowEdge*>::back;
    
};