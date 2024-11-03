#include "Solver.hpp"


class MILPSolver : public Solver{

    IloEnv _env;
    IloModel _model;
    IloCplex _cplex;
    IloArray<IloNumVarArray> _x;
    IloNumVarArray _u;


    void solvemethod(Solution* S);

public:

    MILPSolver(Instance* I);
    ~MILPSolver();

    double gap() {return _cplex.getMIPRelativeGap();};

    Solution* recoversolution();


};