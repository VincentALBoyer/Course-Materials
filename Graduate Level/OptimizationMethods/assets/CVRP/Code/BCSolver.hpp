//Branch and Cut Solver for the TSP with undirected edges
//CPLEX Generic Callbacks implementation (Incumbent)


#pragma once
#include "Solver.hpp"

struct edge{
    int id;
    int i;
    int j;
    double w;
};

class BCSolver : public Solver{

    IloEnv _env;
    IloModel _model;
    IloCplex _cplex;
    IloBoolVarArray _x;

    vector<edge> _edges;

    void solvemethod(Solution* S);

public:

    BCSolver(Instance* I);

    ~BCSolver() { _env.end();};

    double gap() {return _cplex.getMIPRelativeGap();}

    Solution* recoversolution();

   
};

class MyCallBack : public IloCplex::Callback::Function{
    
    Instance* _I;
    IloBoolVarArray _x;
	vector<edge>* _edges;

    void addLazyCuts(const IloCplex::Callback::Context& context);

    void searchHeuristicSolution(const IloCplex::Callback::Context& context);

    void ImproveTour(std::vector<int>& tour, double& length);

	int getEdgeId(int i, int j);

public:
    MyCallBack(Instance* I, IloBoolVarArray& x, vector<edge>* E) 
        : _I(I), _x(x), _edges(E) {}

    // The main callback function
    void invoke(const IloCplex::Callback::Context& context);

};