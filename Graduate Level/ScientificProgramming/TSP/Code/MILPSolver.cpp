#include "MILPSolver.hpp"


MILPSolver::MILPSolver(Instance* I): Solver(I, "MILPSolver"), _env(), _model(_env), _cplex(_model), _x(_env), _u(_env)
{    
    _x = IloArray<IloNumVarArray>(_env, _I->nnodes());
    for(int i=0;i<_I->nnodes();i++){
        _x[i] = IloNumVarArray(_env, _I->nnodes(), 0, 1, IloNumVar::Int);
    }
    _u = IloNumVarArray(_env, _I->nnodes(), 2, _I->nnodes(), IloNumVar::Int);
    
    //Objective function
    IloExpr obj(_env);
    for(int i=0;i<_I->nnodes();i++){
        for(int j=0;j<_I->nnodes();j++){
            if(i!=j){
                obj+=_I->travellingtime(_I->Node(i),_I->Node(j))*_x[i][j];
            }
        }
    }
    _model.add(IloMinimize(_env, obj));
    obj.end();

    //All nodes must be visited
    for(int i=0;i<_I->nnodes();i++){
        IloExpr exprin(_env);
        IloExpr exprout(_env);
        for(int j=0;j<_I->nnodes();j++){
            if(i!=j){
                exprin+=_x[i][j];
                exprout+=_x[j][i];
            }
        }
        _model.add(IloRange(_env, 1, exprin, 1));
        _model.add(IloRange(_env, 1, exprout, 1));
        exprin.end();
        exprout.end();
    }

    //Subtour elimination constraints    
    for(int i=1;i<_I->nnodes();i++){
        for(int j=1;j<_I->nnodes();j++){
            if(i!=j){
                _model.add(_u[i]-_u[j]+1<=(_I->nnodes()-1)*(1-_x[i][j]));
            }
        }
    }

    _cplex.exportModel("Model.lp");
    
}

MILPSolver::~MILPSolver()
{
    _env.end();
}

void MILPSolver::solvemethod(Solution* S)
{
    _cplex.setParam(IloCplex::Param::TimeLimit, _timlim);
    _cplex.setParam(IloCplex::Param::Threads, _threads);
    _cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, _mingap);

    if(S!=NULL){
        //We add MIPStarts
        IloNumVarArray Var(_env); 
        IloNumArray Val(_env);
        for(auto i:*S){
            for(auto j:*S){
                if(i!=j){
                    Var.add(_x[i][j]);
                    Val.add(1);
                }
            }
        } 
        _cplex.addMIPStart(Var, Val);
    }

    _cplex.solve();
    
}

Solution* MILPSolver::recoversolution()
{
    Solution* S = new Solution(_I);
    try{
        int i=0;
        do{
            S->push_back(i);
            for(int j=0;j<_I->nnodes();j++){
                if(i!=j){
                    if(_cplex.getValue(_x[i][j])>0.5){
                        i=j;
                        break;
                    }
                }
            }
        }while(i!=0);
    }
    catch(IloException& e){
        cout << e << endl;
        delete S;
        S = NULL;
    }

    return S;
}