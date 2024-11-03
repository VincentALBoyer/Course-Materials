#include "BCSolver.hpp"


BCSolver::BCSolver(Instance* I): Solver(I, "BCSolver") {

    //We get the undirected edges
    int k = 0;
    for (int i = 0; i < _I->nnodes(); i++) {
        for (int j = 0; j < i; j++) {
            if (i != j) {
                _edges.push_back({ k++, j, i, _I->travellingtime(_I->Node(j), _I->Node(i)) });
            }
        }
    }
    
    _env = IloEnv();
    _model = IloModel(_env);
    _cplex = IloCplex(_model);
    
    _x = IloBoolVarArray(_env, _edges.size());
    for (int i = 0; i < _edges.size(); i++) {
        _x[i] = IloBoolVar(_env);
    }

    //Objective function
    IloExpr obj(_env);
    for (int i = 0; i < _edges.size(); i++) {
        obj += _edges[i].w * _x[i];
    }
    _model.add(IloMinimize(_env, obj));
    obj.end();

    //Constraints
    for (int i = 0; i < _I->nnodes(); i++) {
        IloExpr expr(_env);
        for (int j = 0; j < _edges.size(); j++) {
            if (_edges[j].i == i || _edges[j].j == i) {
                expr += _x[j];
            }
        }
        _model.add(expr == 2);
        expr.end();
    }

}

void BCSolver::solvemethod(Solution* S) {

    if(S!=NULL){
        //We set the initial solution
        IloNumArray Val(_env);
        IloNumVarArray Var(_env);
        for (int i = 0; i < _edges.size(); i++) {
            Var.add(_x[i]);
            Val.add(S->containsUndirectedEdge(_edges[i].i, _edges[i].j));
        }
        _cplex.addMIPStart(Var, Val);
        Val.end();
        Var.end();
    }

    _cplex.setParam(IloCplex::Param::TimeLimit, _timlim);
    _cplex.setParam(IloCplex::Param::Threads, _threads);
    _cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, _mingap);
    

    //We set the callback
    MyCallBack cb(_I, _x, &_edges);
	_cplex.use(&cb, IloCplex::Callback::Context::Id::Candidate | IloCplex::Callback::Context::Id::Relaxation);

    //We solve the model
    _cplex.solve();
}

Solution* BCSolver::recoversolution() {
    cout << "Solver status: " << _cplex.getStatus() << endl;
    Solution* S = new Solution(_I);
    try{
        IloNumArray Val(_env);
        _cplex.getValues(Val, _x);

        //We get the tour starting from node 0
        int i = 0;
        vector<bool> visitededge(_edges.size(), false);
        do {
            //We get the first connected edge
            for (int k = 0; k < _edges.size(); k++) if(!visitededge[k] && Val[k] > 0.5){
            
                if (_edges[k].i == i) {
                    i = _edges[k].j;
                    visitededge[k] = true;
                    break;
                }
                if (_edges[k].j == i) {
                    i = _edges[k].i;
                    visitededge[k] = true;
                    break;
                }
            
            }
            S->push_back(i);
        } while (i != 0);

        Val.end();
    } catch (IloException& e) {
        cout << e << endl;
        delete S;
        S = NULL;
    }
            
    return S;
}


void MyCallBack::addLazyCuts(const IloCplex::Callback::Context& context)
{
    if (!context.inCandidate()) return;

    IloEnv env = context.getEnv();

    IloNumArray Val(env);
    context.getCandidatePoint(_x, Val);

	//We use a rank-based approach to get the subtours
	vector<int> rank(_I->nnodes(), -1);
	
	//We init the ranks
    for (int i = 0; i < _I->nnodes(); i++) {
        rank[i] = -(i+1);
    }

	//We propagate the ranks
    
    auto it = find_if(rank.begin(), rank.end(), [](int r) { return r < 0; });
    while (it!=rank.end()) {
		int i = std::distance(rank.begin(), it);
		int ri = fabs(rank[i]);
		rank[i] = ri;
		for (int j = 0; j < _I->nnodes(); j++) if (i != j) {		 
			int k = getEdgeId(i, j);
			if (Val[k] > 0.5) {
				int rj = fabs(rank[j]);
				if (rj > ri) {
					rank[j] = -ri;
				}
			}	
		}
        it = find_if(rank.begin(), rank.end(), [](int r) { return r < 0; });
	}
   
	//We get the unique ranks
	unordered_set<int> uniqueRanks;
	for (auto r : rank) uniqueRanks.insert(r);

	if (uniqueRanks.size() == 1) return; //No subtour found
   
	//We get the subtours for each rank
    IloRangeArray Cuts(env);
	for (auto r : uniqueRanks) {
		
		//We add a subtour elimination constraint
		IloExpr expr(env);
        int size = 0;
		for (int i = 0; i < _I->nnodes(); i++) if(rank[i]==r){
            size++;
			for (int j = i + 1; j < _I->nnodes(); j++) if (rank[j] == r) {
				int k = getEdgeId(i,j);
				expr += _x[k];
			}
		}
		Cuts.add(expr <= size - 1);
		expr.end();
		
	}

    //We add the constraints
    if (Cuts.getSize() > 0) {
        cout << "Adding " << Cuts.getSize() << " subtour elimination constraints" << endl;
        context.rejectCandidate(Cuts);
    }
    Cuts.end();

    Val.end();
}

void MyCallBack::searchHeuristicSolution(const IloCplex::Callback::Context& context)
{
	if (!context.inCandidate() && !context.inRelaxation()) return;

    if (context.inRelaxation()) {
        //long Node_ID = context.getLongInfo(IloCplex::Callback::Context::Info::NodeUID);
        CPXLONG Node_LEFT = context.getLongInfo(IloCplex::Callback::Context::Info::NodesLeft);
        if (Node_LEFT > 1 && RAND01() > 0.1) return;
    }
    else return;

	IloEnv env = context.getEnv();
    IloNumArray Var(env);
	if (context.inRelaxation()) context.getRelaxationPoint(_x, Var);
	else context.getCandidatePoint(_x, Var);

	//We greedily build a tour by adding edges with the highest value x connected to both extremity to the tour
	//If the edge is already in the tour we skip it
	//We resolve equalities weights by adding the connected edge with the lowest weight
	vector<int> tour;
	vector<bool> visitednode(_I->nnodes(), false);
	int n = _I->nnodes();

	//We get the first edge
	auto it = max_element(_edges->begin(), _edges->end(), [&Var](const edge& e1, const edge& e2) { return Var[e1.id] < Var[e2.id] || (Var[e1.id] == Var[e2.id] && e1.w < e2.w); });
	tour.push_back(it->i);
	tour.push_back(it->j);
	n -= 2;

    while (n != 0)
    {
        //we get the extremities of the tour
        int i = tour.front();
        int j = tour.back();
        visitednode[i] = true;
        visitednode[j] = true;

        //We get the best edge connected to one of the extremities
        int bestedgeid = -1;
        double bestedgeval = -1;
        double bestedgew = -1;
		for (int l = 0; l < _I->nnodes(); l++) if (!visitednode[l]){
			int ki = getEdgeId(i, l);
			double wi = (*_edges)[ki].w;
		    if (Var[ki] > bestedgeval || (Var[ki] == bestedgeval && wi < bestedgew)) {
			    bestedgeid = ki;
			    bestedgeval = Var[ki];
			    bestedgew = wi;
		    }

			int kj = getEdgeId(j, l);
			double wj = (*_edges)[kj].w;
			if (Var[kj] > bestedgeval || (Var[kj] == bestedgeval && wj < bestedgew)) {
				bestedgeid = kj;
				bestedgeval = Var[kj];
				bestedgew = wj;
			}
		}

        if (bestedgeid == -1) 
            throw exception("No edge found (Heuristic)");

        //We add the edge to the tour
        if ((*_edges)[bestedgeid].i == i) tour.insert(tour.begin(), (*_edges)[bestedgeid].j);
        else if ((*_edges)[bestedgeid].j == i) tour.insert(tour.begin(), (*_edges)[bestedgeid].i);
        else if ((*_edges)[bestedgeid].i == j) tour.push_back((*_edges)[bestedgeid].j);
        else tour.push_back((*_edges)[bestedgeid].i);
        n--;
    }

	//We compute the length of the tour
	double length = 0;
    for (int i = 0; i < tour.size(); i++) {
        int j = (i + 1) % tour.size();
		int k = getEdgeId(tour[i], tour[j]);
		length += (*_edges)[k].w;
    }

    //We try to improve the lenght by swapping contiguous nodes in the tour
    //ImproveTour(tour, length);

	//We add the solution to the incumbent
    double incumbent = context.getIncumbentObjective();
	if (incumbent > length) {
		vector<bool> isEdgeInTour(_edges->size(), false);
		for (int i = 0; i < tour.size(); i++) {
			int j = (i + 1) % tour.size();
			auto it = find_if(_edges->begin(), _edges->end(), [&tour, i, j](const edge& e) { 
                return (e.i == tour[i] && e.j == tour[j]) || (e.j == tour[i] && e.i == tour[j]); 
                });
			if (it != _edges->end()) isEdgeInTour[it->id] = true;
			else throw exception("Edge not found (Heuristic)");
		}

		IloNumVarArray X(env);
		IloNumArray XVal(env);
		for (int i = 0; i < _edges->size(); i++) {
			X.add(_x[i]);
			XVal.add(isEdgeInTour[i]);
		}
		context.postHeuristicSolution(X, XVal, length, IloCplex::Callback::Context::SolutionStrategy::CheckFeasible);
		X.end();
		XVal.end();
    }
	

    


    Var.end();
}

void MyCallBack::ImproveTour(std::vector<int>& tour, double& length)
{
    bool improved = true;
    while (improved) {
        improved = false;
        for (int i = 0; i < tour.size(); i++) {
            int j = (i + 1) % tour.size();
            int k = (i + 2) % tour.size();
            int l = (i + 3) % tour.size();
            double delta = 0;
            delta -= (*_edges)[getEdgeId(tour[i], tour[j])].w;
            //delta -= (*_edges)[getEdgeId(tour[j], tour[k])].w;
            delta -= (*_edges)[getEdgeId(tour[k], tour[l])].w;
            delta += (*_edges)[getEdgeId(tour[i], tour[k])].w;
            delta += (*_edges)[getEdgeId(tour[j], tour[l])].w;
            if (delta < 0) {
                swap(tour[j], tour[k]);
                length += delta;
                improved = true;
            }
        }
    }
}

int MyCallBack::getEdgeId(int i, int j){
	int ii = max(i, j);
	int jj = min(i, j);
	return (ii * (ii - 1)) / 2 + jj;
}

void MyCallBack::invoke(const IloCplex::Callback::Context& context) {

	addLazyCuts(context);

	searchHeuristicSolution(context);
}