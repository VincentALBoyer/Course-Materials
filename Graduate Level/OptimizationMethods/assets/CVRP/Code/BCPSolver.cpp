#include "BCPSolver.hpp"



unsigned int BBNode::_index = 0;




void CGMasterProblem::setBounds(BBNode* N) {
	for (auto x : _X) if(x->getType() == CGEnv::EdgeVar ) {
		x->setLB(0.0);
		x->setUB(1.0);
	}

	//static vector<int> _lastIndexMod;
	//for (auto k : _lastIndexMod) if (k >= 0) {
	//	_X[k]->setLB(0.0);
	//	_X[k]->setUB(1.0);
	//}
	//_lastIndexMod.clear();

	BBNode* node = N;
	while (node != NULL) {
		//if (node->getVar() != NULL) _lastIndexMod.push_back(node->getVar()->getId());
		node->setVarBound();
		//cout << node->getVar()->getName() << "(" << node->getVar()->getLB() << "," << node->getVar()->getUB() << ")" << endl;
		node = node->getParent();
	}

	vector<CGNumVar*> prunedVars = pruneNode(N);
	//for (auto X : prunedVars) {
	//	_lastIndexMod.push_back(X->getId());
	//}


	


}

vector<CGNumVar*> CGMasterProblem::pruneNode(BBNode* N)
{
	//We retreive all the Up BBEdges from the root to the node
	vector<CGEdge> edges;
	BBNode* node = N;
	while (node->getParent() != NULL) {
		if (node->getNodeDir() == CGEnv::Up) {
			CGNumVar* x = node->getVar();
			edges.push_back({ x->geti(), x->getj(), 0 });
		}
		node = node->getParent();
	}

	//We reconstruct all partial tours from the set of edges
	vector< vector<int> > tours;
	while (!edges.empty()) {
		vector<int> tour;
		tour.push_back(edges.back().i);
		tour.push_back(edges.back().j);
		edges.pop_back();
		bool hasChanged = true;
		while (hasChanged) {
			hasChanged = false;
			for (auto it = edges.begin(); it != edges.end(); it++) {
				if (tour.back() == it->i) {
					tour.push_back(it->j);
					edges.erase(it);
					hasChanged = true;
					break;
				}
				else if (tour.back() == it->j) {
					tour.push_back(it->i);
					edges.erase(it);
					hasChanged = true;
					break;
				}
				else if (tour.front() == it->i) {
					tour.insert(tour.begin(), it->j);
					edges.erase(it);
					hasChanged = true;
					break;
				}
				else if (tour.front() == it->j) {
					tour.insert(tour.begin(), it->i);
					edges.erase(it);
					hasChanged = true;
					break;
				}
			}
		}
		tours.push_back(tour);
	}

	//We forbids the edges that will yield a subtour
	int n = _env.getInstance()->nnodes();
	vector<CGNumVar*> forbidden;
	for (auto& t : tours) if(t.size()>2 && t.size()<n ) {
		int i = t.back();
		int j = t.front();
		//We forbid the edge
		CGNumColumn* col = getCol(i, j);
		if (col->isInModel()) {
			CGNumVar* x = col->getVar();
			x->setLB(0.0);
			x->setUB(0.0);
			forbidden.push_back(x);
		}
		
	}

	return forbidden;

}

CGMasterProblem::CGMasterProblem(Instance* I) : _env(I)
{
	_maxColPerRound = 100;

	//We generate all the columns
	_ncols = I->nnodes() * (I->nnodes() -1) / 2;
	_columns = new CGNumColumn * [_ncols];
	int index = 0;
	for (int i = 0; i < I->nnodes(); i++) {
		for (int j = 0; j < i; j++) {
			_columns[getColIndex(i, j)] = new CGNumColumn(i, j, getObjCoeff(i, j));
			assert(index >= 0 && index < _ncols);
			assert(index == getColIndex(i, j));
			index++;
			//cout << "Column " << getColIndex(i, j) << ": " << i << "-" << j << " : " << getObjCoeff(i, j) << endl;
		}
	}
	assert(index == _ncols);

	

	//for (int c = 0; c < _ncols; c++) 
	//	cout << _columns[c]->geti() << "-" << _columns[c]->getj() << " : " << _columns[c]->getObj() << endl;
	
	_obj = IloObjective(_env, 0, IloObjective::Minimize);

	_model = IloModel(_env);
	_cplex = IloCplex(_model);

	_model.add(_obj);

	//We add the degree constraints
	for (int i = 0; i < I->nnodes(); i++) {
		CGConstraintNodeDegree* constraint = new CGConstraintNodeDegree(_env, i);
		addConstraint(constraint);
	}

	//CGConstraint* cut = new CGConstraintMaxNEdges(_env, I->nnodes());
	//addConstraint(cut);

	//_cplex.exportModel("M.lp");
	//We turn off the output
	_cplex.setOut(_env.getNullStream());
	
	//_cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Algorithm::Primal);



	

}

CGMasterProblem::~CGMasterProblem()
{
	for (auto x : _X) delete x;
	for (auto c : _constraints) delete c;
	for (int i = 0; i < _ncols; i++) delete _columns[i];
	delete[] _columns;
}

void CGMasterProblem::addVar(int i, int j)
{
	CGNumColumn* col = getCol(i, j);
	addCol(col);
}

void CGMasterProblem::addConstraint(CGConstraint* constraint)
{
	//We try to merge the constraint with the existing ones: Bag Idea
	//for (auto R : _constraints) if (R->merge(constraint)) {
	//	//We update all the columns
	//	for (int c = 0; c < _ncols; c++) {
	//		CGNumColumn* col = _columns[c];
	//		double coef = R->getCoefficient(col->geti(), col->getj());
	//		col->add(R->getId(), coef);
	//		if (col->isInModel()) R-> setLinearCoef(*(col->getVar()), coef);
	//	}
	//	delete constraint;
	//	return;
	//}



	CGEnv::ConstraintType type = constraint->getType();

	//We add the slack variable if needed
	double lb = constraint->getLB();
	double ub = constraint->getUB();
	if (fabs(ub - lb) > CGEnv::EPSILON) {
		if (lb > -IloInfinity && ub < IloInfinity)
			throw runtime_error("Invalid constraint bounds");
		if (lb > -IloInfinity) {
			CGNumVar* slack =  new CGNumVar(_env);
			constraint->setLinearCoef(*slack, -1.0);
			_model.add(*slack);
			_X.push_back(slack);
			constraint->setLB(lb);
			constraint->setUB(lb);

		}
		else if (ub < IloInfinity) {
			CGNumVar* slack = new CGNumVar(_env);
			constraint->setLinearCoef(*slack, 1.0);
			_model.add(*slack);
			_X.push_back(slack);
			constraint->setLB(ub);
			constraint->setUB(ub);
		}
	}


	

	//assert(type != CGEnv::NotSpecified);
	constraint->setId(unsigned int(_constraints.size()));
	_constraints.push_back(constraint);
	_model.add(*constraint);
	int k = int(_constraints.size() - 1);

	//We only update columns not in model (Be carefull for debugging)
	//The other columns must have been updated when the constraint was created
	for (int c = 0; c < _ncols; c++) {
		CGNumColumn* col = _columns[c];
		//if (col->isInModel()) continue;
		double coef = 0.0;
		if (constraint->getType() == CGEnv::Gomory) {
			CGConstraintGomory* gomory = dynamic_cast<CGConstraintGomory*>(constraint);
			coef = gomory->getCoefficient( col);
		}
		else coef = constraint->getCoefficient(col->geti(), col->getj());
		//if (col->isInModel() && constraint->getType() == CGEnv::Gomory)
		//	coef = constraint->retrieveCoefficient(col->getVar());
		col->add(k, coef);

	}

	constraint->clean();

	//Instance* I = _env.getInstance();
	//for (int i = 0; i < I->nnodes(); i++)
	//	for (int j = 0; j < i; j++) getCol(i, j);

	

}

void CGMasterProblem::addConstraints(vector<CGConstraint*>& constraints)
{
	for (auto c : constraints) addConstraint(c);
}

CGEnv::SolverStatus CGMasterProblem::solve(BBNode* N, double UB)
{
	static int offset = 0;

	offset = offset % _ncols;

	Instance* I = _env.getInstance();

	if (N == NULL)
		throw runtime_error("Null Node passed to master problem");

	//UB = 51200;
		
	//We reset the bounds
	setBounds(N);
	
	
	//Columns + Shuffle
	//vector<CGNumColumn*> C;
	//for (int c = 0; c < _ncols; c++) if(!_columns[c]->isInModel()) C.push_back(_columns[c]);
	//sort(C.begin(), C.end(), [&](CGNumColumn* a, CGNumColumn* b) {
	//	return a->getObj() < b->getObj();
	//	});

	//We solve the master problem with column generation
	auto start = chrono::high_resolution_clock::now();
	bool hasNewColumns = true;
	int nvarsadded = 0;
	double targetObj = N->getObjVal();
	double currentObj = DBL_MAX;
	while (hasNewColumns) {

		//_cplex.exportModel("M.lp");
		_cplex.solve();

		//ºWe check the status
		if (_cplex.getStatus() != IloAlgorithm::Optimal) return CGEnv::SolverStatus::Infeasible;
		currentObj = _cplex.getObjValue();

		if (fabs(currentObj- targetObj)<CGEnv::EPSILON ) 
			break;

		//We get the duals
		IloNumArray duals = getDuals();
		double sum = -DBL_MAX;
		for (int k = 0; k < _ncols; k++) {
			CGNumColumn* col = _columns[k];
			if(!col->isInModel())
				sum = fmax(sum, duals[col->geti()] + duals[col->getj()]);
		}
		for (int k = I->nnodes(); k < duals.getSize(); k++) if(duals[k]>CGEnv::EPSILON) sum += duals[k];


		vector<CGNumColumn*> candidates;
		for (int k = 0; k < _ncols; k++, offset++) {
			CGNumColumn* col = _columns[(offset)%_ncols];
			if (col->getObj() - sum > -CGEnv::EPSILON || col->isInModel()) continue;
			double rc = col->getReducedCost(duals);
			if (rc < -CGEnv::EPSILON) {
				candidates.push_back(col);
				if (candidates.size() >= 100) break;
			}
		}

		//We add the new variables
		for (auto& e : candidates) {
			addCol(e);
		}

		hasNewColumns = !candidates.empty();

		nvarsadded += int(candidates.size());

		duals.end();

	}

	N->setObjVal(currentObj);

	auto end = chrono::high_resolution_clock::now();

	//We Display the log
	//cout << "MP: Obj=" <<_cplex.getObjValue() << ", #addedvars = " << nvarsadded << ", #vars = " << _X.size() << ", #csts = " << _constraints.size() 
	//	<< ", Time = " << chrono::duration_cast<chrono::seconds>(end - start).count() << "s" << endl;

	//We pruned variables locally
	IloNumArray duals = getDuals();
	for (int k = 0; k < _ncols; k++) {
		CGNumColumn* col = _columns[k];
		if (!col->isInModel() || col->getVar()->getUB()>0.5) continue;
		double rc = _cplex.getReducedCost(*col->getVar());
		if (rc >= UB - currentObj) {
			N->setAsPruned(col->getVar());
		}
	}
	duals.end();
	

	return CGEnv::Optimal;
}

IloNumArray CGMasterProblem::getVarVals()
{
	// TODO: insert return statement here
	IloNumVarArray vars(_env);
	IloNumArray vals(_env);

	for (auto x : _X) vars.add(*x);

	_cplex.getValues(vals, vars);
	vars.end();
	return vals;

}

IloNumArray CGMasterProblem::getDuals()
{
	// TODO: insert return statement here
	IloRangeArray ranges(_env);
	IloNumArray duals(_env);
	for (int k = 0; k < _constraints.size(); k++) ranges.add(*_constraints[k]);
	_cplex.getDuals(duals, ranges);
	ranges.end();

	return duals;
}


vector<double> CGMasterProblem::getInvRow(int k)
{
	// Get the CPLEX environment and problem pointer
	CPXENVptr envPtr = _cplex.getImpl()->getCplexEnv();
	CPXLPptr lpPtr = _cplex.getImpl()->getCplexLp();

	// Specify the row index you want to retrieve
	std::vector<double> row(_cplex.getNcols());

	// Call CPXbinvrow to get the i-th row of the basis inverse
	int status = CPXbinvarow(envPtr, lpPtr, k, row.data());
	assert(status == 0);

	return row;
}

vector<double> CGMasterProblem::getInvBasisRow(int k)
{
	// Get the CPLEX environment and problem pointer
	CPXENVptr envPtr = _cplex.getImpl()->getCplexEnv();
	CPXLPptr lpPtr = _cplex.getImpl()->getCplexLp();

	int numRows = _cplex.getNrows();
	std::vector<double> row(numRows);
	std::vector<double> rhs(numRows);
	if (numRows != _constraints.size())
		throw runtime_error("Invalid number of rows");

	// Call CPXbinvrow to get the i-th row of the basis inverse
	int status = CPXXbinvrow(envPtr, lpPtr, k, row.data());
	assert(status == 0);
	CPXgetrhs(envPtr, lpPtr, rhs.data(), 0, numRows - 1);

	return row;
}

double CGMasterProblem::getInvRHS(int k)
{
	// Get the CPLEX environment and problem pointer
	CPXENVptr envPtr = _cplex.getImpl()->getCplexEnv();
	CPXLPptr lpPtr = _cplex.getImpl()->getCplexLp();

	int numRows = _cplex.getNrows();
	std::vector<double> row(numRows);
	std::vector<double> rhs(numRows);
	if (numRows != _constraints.size())
		throw runtime_error("Invalid number of rows");

	// Call CPXbinvrow to get the i-th row of the basis inverse
	int status = CPXXbinvrow(envPtr, lpPtr, k, row.data());
	assert(status == 0);
	CPXgetrhs(envPtr, lpPtr, rhs.data(), 0, numRows-1);

	//We do the dot product
	double invrhs = 0.0;
	for (int i = 0; i < numRows; i++) {
		invrhs += row[i] * rhs[i];
	}

	return invrhs;
}

pair<int*, double*> CGMasterProblem::getBHead()
{
	// Get the CPLEX environment and problem pointer
	CPXENVptr envPtr = _cplex.getImpl()->getCplexEnv();
	CPXLPptr lpPtr = _cplex.getImpl()->getCplexLp();

	int* head = new int[_cplex.getNrows()];
	double* x = new double[_cplex.getNrows()];

	CPXgetbhead(envPtr, lpPtr, head, x);

	return make_pair(head, x);
}

IloCplex::BasisStatusArray CGMasterProblem::getStatus()
{
	IloCplex::BasisStatusArray cstat(_env);
	IloNumVarArray Var(_env);
	for (auto x : _X) Var.add(*x);
	
	_cplex.getBasisStatuses(cstat, Var);
	Var.end();
	return cstat;
}

double CGMasterProblem::getReducedCost(int i, int j, IloNumArray& duals) {
	
	CGNumColumn* C = getCol(i, j);
	return getReducedCost(C, duals);

}

double CGMasterProblem::getReducedCost(CGNumColumn* col, IloNumArray& duals)
{
	CGNumColumn& C = *col;
	double rc = C.getObj();
	for (auto& it : C) rc -= duals[it.first] * it.second;

	return rc;
}

void CGMasterProblem::addCol(CGNumColumn* col)
{
	int nextid = int(_X.size());
	IloNumColumn iloCol(_env);
	iloCol += _obj(col->getObj());
	//for (int k = 0; k < _constraints.size(); k++) if(fabs(C[k])>CGEnv::EPSILON) iloCol += (*_constraints[k])(C[k]);
	for (auto& it : *col) iloCol += (*_constraints[it.first])(it.second);
	CGNumVar* x = new CGNumVar(_env, iloCol, col->geti(), col->getj());
	assert(nextid == x->getId());
	_X.push_back(x);
	_model.add(*x);

	col->setInModel(x);
	col->clear();

	iloCol.end();
}

double CGMasterProblem::getDistanceToCut(CGConstraint& R, IloNumArray& vals)
{
	return 0.5;
	IloExpr expr(R.getExpr());
	IloNum lhs = R.getLB();
	IloNum rhs = R.getUB();

	IloNum val = _cplex.getValue(expr);

	double diff = fmax((lhs - val), (val - rhs));
	assert(diff >= 0.0);
	expr.end();
	return diff;

	IloExpr::LinearIterator it = expr.getLinearIterator();
	double qsum = 0.0;
	while (it.ok()) {
		if (it.getVar().getUB() > CGEnv::EPSILON)
			qsum += pow(it.getCoef(), 2);
		++it;
	}
	if (qsum <  CGEnv::EPSILON) {
		cout << R << endl;
	}
	if (qsum < CGEnv::EPSILON) return NULL;
	assert(qsum > CGEnv::EPSILON);
	qsum = sqrt(qsum);

	double dist = diff / qsum;
	//cout << _R->getName() <<  ": " << diff << " - " << qsum << " - " << dist << endl;
	//assert(fabs(diff) >= fabs(dist));

	expr.end();
	return dist;
}

void BCPSolver::solvemethod(Solution* S)
{


	//We initialise the cout header
	int w = 15;
	cout << setw(w) << "Left" << setw(w) << "NodeId" << setw(w) << "NodeO" << setw(w) << "BObjVal" << setw(w) << "BBound" << setw(w) << "Gap" << endl;

	//We initialize the Master Problem
	_bestobjval = DBL_MAX;
	if (_master->getNumVars() == 0) initMasterProblem();
	//addAllVars();

	//We create the root node
	BBNode* root = new BBNode(_master->getEnv(),0);
	_BBNodes.push_back(root);

	//We Enter the cutting Phase
	_cuttingplanesolver->solve(root);

	//We start the branch and bound
	unsigned int iter = 0;
	while (!_BBNodes.empty() && gap()>CGEnv::EPSILON && timeLeft()>CGEnv::EPSILON) {
		iter++;
		//We get the best node
		BBNode* node = _BBNodes.back();
		_BBNodes.pop_back();
		if(!_BBNodes.empty() && node->getObjVal()>_BBNodes.back()->getObjVal())
			throw runtime_error("Invalid node order");

		//if (node->getDepth() < 2) _cuttingplanesolver->solve(node);

		//We get the current best bound
		double bestbound = node->getObjVal();

		//We solve the master problem
		CGEnv::SolverStatus status = _master->solve(node, _bestobjval );

		if (status == CGEnv::Infeasible && node->getDepth() == 0)
			throw runtime_error("Infeasible problem at the root");
		else if (status == CGEnv::Infeasible || _master->getObjVal() >= _bestobjval) {
			node->setStatus(CGEnv::InfeasibleNode);
		}
		else {
			//node->setObjVal(_master->getObjVal());
			node->setStatus(CGEnv::FeasibleNode);
		}

		if (node->getStatus() == CGEnv::FeasibleNode) {
			//We get the solution
			IloNumArray vals = _master->getVarVals();

			//A Quick Check if the solution is consistent with node varstatus
			if (node->getDepth() > 0) {
				double valk = vals[node->getVar()->getId()];
				if (node->getVar()->getName() != _master->getVar(node->getVar()->getId())->getName()) {
					cout << node->getVar()->getName() << "!=" << _master->getVar(node->getVar()->getId())->getName() << endl;
					throw runtime_error("Inconsistent variable name: " + string(node->getVar()->getName()) + "!=" + string(_master->getVar(node->getVar()->getId())->getName()));
				}
				if (node->getNodeDir() == CGEnv::Down && valk > CGEnv::EPSILON)
					throw runtime_error("Inconsistent solution (D): " + to_string(valk));
				else if (node->getNodeDir() == CGEnv::Up && valk < 1 - CGEnv::EPSILON)
					throw runtime_error("Inconsistent solution (U): " + to_string(valk));
				//else throw runtime_error("Unknown node varstatus");
			}

			//We get the fractional variables
			vector<CGNumVar*> fracvarsindex;
			for (int i = 0; i < _master->getNumVars(); i++) 
				if (CGEnv::isFractional(vals[i]) && _master->getVar(i)->getType()==CGEnv::EdgeVar)
					fracvarsindex.push_back(_master->getVar(i));

			//We check if the solution is integer feasible
			if (fracvarsindex.empty()) {
				vector<CGConstraint*> Cuts;
				//_cuttingplanesolver->getInfeasiblePathCuts(vals, Cuts, _bestobjval);
				if (!_cuttingplanesolver->getSECCuts(vals, Cuts)) {
					node->setStatus(CGEnv::IntegerNode);
					cout << "Integer Node: " << node->getObjVal() << endl;
				}
				else {
					node->setStatus(CGEnv::IntegerInfeasibleNode);
					_master->addConstraints(Cuts);
				}
			}
			else if (fracvarsindex.size() < fmax(10, 0.05 * _master->getNumConstraints())) {
				applyHeuristic(vals);
			}

			if (node->getStatus() == CGEnv::FeasibleNode) {
			
				//Custom branching
				CGNumVar* X = NULL;
				double maxfrac = -1;
				//auto candidates = getBranchingCandidates(vals);
				for (auto x : fracvarsindex) {
					double val = vals[x->getId()];
					double frac = fmin(val, 1.0 - val) * _I->travellingtime(_I->Node(x->geti()), _I->Node(x->getj()));
					if (frac > maxfrac) {
						maxfrac = frac;
						X = x;
					}
					if (maxfrac < CGEnv::EPSILON) break;
				}

				if (X == NULL) throw runtime_error("No fractional variable found");
				//if(!CGEnv::isFractional(vals[X->getId()]))
				//	throw runtime_error("Invalid variable status: " + to_string(vals[X->getId()]));
				//if (maxfrac < CGEnv::EPSILON)
				//	throw runtime_error("Invalid maxfrac");
				//We create the left and right nodes
				node->generateChildren(X);

				//We insert the children in the list
				addNode(node->getLeft());
				addNode(node->getRight());
			}
			else if (node->getStatus() == CGEnv::IntegerInfeasibleNode) {
				//We reinsert the node in the list
				addNode(node);
			}
			else if (node->getStatus() == CGEnv::IntegerNode) {
				updateBestTour(vals);  //Not Covered by heuristic (Bug?)
			}
			else throw runtime_error("Unknown node status");

			vals.end();
		}

		
		
		

		//We print some logs: gap, bestobjval, bestbound
		if (iter ==1 || iter % 100 == 0) {
			string Info = to_string(node->getId());
			if (node->getStatus() == CGEnv::FeasibleNode) Info += "";
			else if (node->getStatus() == CGEnv::InfeasibleNode) Info = "Inf";
			else if (node->getStatus() == CGEnv::IntegerNode) Info += "*";
			else if (node->getStatus() == CGEnv::IntegerInfeasibleNode) Info += "IInf";
			else if (node->getStatus() == CGEnv::PendingNode) Info += "";
			else throw runtime_error("Unknown node status");
			cout << setw(w) << _BBNodes.size() << setw(w) << Info << setw(w) << node->getObjVal() << setw(w) << _bestobjval << setw(w) << bestbound << setw(w) << gap() << endl;
			
		}
		

		//We clean the infeasible node
		if (node->getStatus() == CGEnv::InfeasibleNode || node->getStatus() == CGEnv::IntegerNode) 
			delete node;
	}

	cout << "Best Obj Val: " << _bestobjval << endl;
	//cout << "Best Tour: ";
	//for (int i = 0; i < _besttour.size(); i++) cout << _besttour[i] << " ";
	//cout << endl;


	
}

void BCPSolver::initMasterProblem(const int maxGenTours, const int k)
{
	//We generate some random path with nearest neighbor
	//And add the unique edges to the master problem
	set< CGEdge > edges;
	for (int i = 0; i < maxGenTours; i++) {

		//We generate a random tour with k-nearest neighbor
		vector<int> tour;
		vector<bool> visited(_I->nnodes(), false);
		double* distances = new double[_I->nnodes()];
		int current = rand() % _I->nnodes();
		tour.push_back(current);
		visited[current] = true;
		for (int j = 1; j < _I->nnodes(); j++) {
			vector<int> nearest;
			
			for (int l = 0; l < _I->nnodes(); l++)  if (!visited[l]) {
				nearest.push_back(l);
				distances[l] = _I->travellingtime(_I->Node(tour.back()), _I->Node(l));
			}

			//We get the k-nearest neighbors 
			vector<int> knearest;
			for (int n = 0; n < min(k, (int)nearest.size()) ; n++) {
				auto it = min_element(nearest.begin(), nearest.end(), [&](int i, int j) {
					return distances[i] < distances[j];
					});
				knearest.push_back(*it);
				nearest.erase(it);
			}

			//We randomly select a nearest neighbor
			int next = knearest[rand() % knearest.size()];
			tour.push_back(next);
			visited[next] = true;

		}
		delete[] distances;


		//We collect the unique edges
		for (int j = 0; j < _I->nnodes(); j++) {
			int o = max(tour[j], tour[(j + 1) % _I->nnodes()]);
			int d = min(tour[j], tour[(j + 1) % _I->nnodes()]);
			edges.insert({ d, o, 0 });
		}

		//We compute the length of the tour
		addFeasibleTour(tour);


	}

	//We add the unique edges to the master problem
	for (auto& e : edges) {
		_master->addVar(e.i, e.j);
	}
}

void BCPSolver::initMasterProblemRandom(const int maxGenTours)
{
	//We generate some random path with nearest neighbor
	//And add the unique edges to the master problem
	set< CGEdge > edges;
	//vector<int> tour;
	//for (int i = 0; i < _I->nnodes(); i++) tour.push_back(i);
	//for (int i = 0; i < maxGenTours; i++) {
	//	//We shuffle the tour
	//	random_shuffle(tour.begin(), tour.end());

	//	//We collect the unique edges
	//	for (int j = 0; j < _I->nnodes(); j++) {
	//		int o = max(tour[j], tour[(j + 1) % _I->nnodes()]);
	//		int d = min(tour[j], tour[(j + 1) % _I->nnodes()]);
	//		edges.insert({ d, o, 0 });
	//	}
	//}

	//We only add the two shrtest edges for each node
	for (int i = 0; i < _I->nnodes(); i++) {
		vector<int> nearest;
		for (int j = 0; j < _I->nnodes(); j++) if (i != j) nearest.push_back(j);
		sort(nearest.begin(), nearest.end(), [&](int i, int j) {
			return _I->travellingtime(_I->Node(i), _I->Node(j)) < _I->travellingtime(_I->Node(j), _I->Node(i));
			});
		for (int k = 0; k < 2; k++) {
			int j = nearest[k];
			int o = max(i, j);
			int d = min(i, j);
			edges.insert({ d, o, 0 });
		}
	}

	//We add the unique edges to the master problem
	for (auto& e : edges) {
		_master->addVar(e.i, e.j);
	}
		
}

void BCPSolver::addAllVars()
{
	//For debugging purposes
	for (int i = 0; i < _I->nnodes(); i++) {
		for (int j = 0; j < i; j++) {
			_master->addVar(i, j);
		}
	}
}

bool BCPSolver::updateBestTour(IloNumArray& vals)
{
	if (!CGEnv::isInteger(vals)) return false;
	//We get the current tour starting with node 0
	vector<int> tour;
	vector<bool> visited(_I->nnodes(), false);
	int current = 0;

	tour.push_back(current);
	visited[current] = true;
	bool nextfound = true;
	while(nextfound && tour.size()<_I->nnodes()) {
		nextfound = false;
		for (int k = 0; k < vals.getSize() && !nextfound; ++k) {
			if (vals[k] > 0.5) {
				CGNumVar* x = _master->getVar(k);
				if (x->geti() == current && !visited[x->getj()]) {
					current = x->getj();
					nextfound = true;
				}
				if (x->getj() == current && !visited[x->geti()]) {
					current = x->geti();
					nextfound = true;
				}
			}
		}
		if (nextfound) {
			tour.push_back(current);
			visited[current] = true;
		}
	}

	if (tour.size() != _I->nnodes()) 
		throw runtime_error("Invalid tour size: " + to_string(tour.size()));

	//CGConstraint* Cut = new CGConstraintInfPath(_master->getEnv(), tour);
	//_master->setupConstraint(Cut);
	//_master->addConstraint(Cut);

	return addFeasibleTour(tour);
	
}

bool BBCuttingPlaneSolver::getSECCuts(IloNumArray& vals, vector<CGConstraint*>& Cuts)
{
	if (!CGEnv::isInteger(vals)) return false;
	
	//We use a rank-based approach to get the subtours
	vector<int> rank(_I->nnodes(), -1);

	//We init the ranks
	for (int i = 0; i < _I->nnodes(); i++) {
		rank[i] = -(i + 1);
	}

	size_t initSize = Cuts.size();

	//We propagate the ranks
	auto it = find_if(rank.begin(), rank.end(), [](int r) { return r < 0; });
	while (it != rank.end()) {
		int i = std::distance(rank.begin(), it);
		assert(rank[i] < 0);
		int ri = fabs(rank[i]);
		rank[i] = ri;
		for (int k = 0; k < vals.getSize(); ++k) if (vals[k] > 0.5) {
			CGNumVar* x = _master->getVar(k);
			int xi = x->geti();
			int xj = x->getj();
			if (xi == i || xj == i) {
				int j = xi == i ? xj : xi;
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

	if (uniqueRanks.size() == 1) {
		return false; //No subtour found
		//updateBestTour(vals);
	}

	//We get the subtours for each rank
	for (auto r : uniqueRanks) {

		//We add a subtour elimination constraint
		vector<int> subtour;
		for (int i = 0; i < _I->nnodes(); i++) if (rank[i] == r) subtour.push_back(i);
		CGConstraintSubtourElimination* cut = new CGConstraintSubtourElimination(_master->getEnv(), subtour);
		_master->setupConstraint(cut);
		Cuts.push_back(cut);

	}

	cout << "#SEC Found: " << uniqueRanks.size() << endl;

	return initSize < Cuts.size();

}

bool BBCuttingPlaneSolver::getGomoryCuts(IloNumArray& vals, vector<CGConstraint*>& Cuts)
{

	if (CGEnv::isInteger(vals)) return false;

	//return false;
	//Preprocessing
	vector<CGConstraint*> CSTS(_master->getNumConstraints());
	for (int k = 0; k < _master->getNumConstraints(); ++k) {
		CSTS[k] = _master->getConstraint(k);
	}

	
	size_t initSize = Cuts.size();

	vector< pair<CGConstraint*, double> > CutsCollected;
	//IloNumArray duals = _master->getDuals();

	vector<double> RHS(_master->getNumConstraints());
	vector<int> kindex(_master->getNumConstraints());
	for (int k = 0; k < _master->getNumConstraints(); k++) {
		RHS[k] = _master->getInvRHS(k);
		kindex[k] = k;
	}

	auto bhead = _master->getBHead();
	//std::sort(kindex.begin(), kindex.end(), [&](int i, int j) {return CGEnv::FractionalPart(RHS[i]) > CGEnv::FractionalPart(RHS[j]); });
	//std::sort(kindex.begin(), kindex.end(), [&](int i, int j) {return fabs(0.5 - RHS[i]) < fabs(0.5 - RHS[j]); });

	vector<double> BVarVals(_master->getNumConstraints());
	for (int k = 0; k < _master->getNumConstraints(); k++) {
		int i = bhead.first[k];
		if (i >= 0) {
			CGNumVar* x = _master->getVar(i);
			BVarVals[k] = _master->getVarVal(*x);
		}
		else BVarVals[k] = 0.0;
	}
	std::sort(kindex.begin(), kindex.end(), [&](int i, int j) {return fabs(0.5 - BVarVals[i]) < fabs(0.5 - BVarVals[j]); });

	
	//for (int k = 0; k < _master->getNumConstraints(); k++) {
	//	int i = bhead.first[k];
	//	if (i >= 0) {
	//		CGNumVar* x = _master->getVar(i);
	//		cout << x->getName() << " " << bhead.second[k] << " " << RHS[k] << " " << _master->getVarVal(*x) << endl;
	//	}
	//	
	//}

	//We generate the Gomory cuts for each constraint
	for (int k : kindex) if(CutsCollected.size()<1){

		if (bhead.first[k] < 0) continue;
		CGNumVar* BasicVar = _master->getVar(bhead.first[k]);
		double BasicVal = bhead.second[k];

		if (CGEnv::isInteger(BasicVal) || BasicVar->getType()==CGEnv::VarType::SlackVar ) continue;

		//if (_master->getConstraint(k)->getType() == CGEnv::Gomory) continue;

		CGConstraint* R = _master->getConstraint(k);
		
		if (R->getType() != CGEnv::NotSpecified) {

			double rhs = RHS[k];// _master->getInvRHS(k);
			vector<double> row = _master->getInvRow(k);
			auto cstat = _master->getStatus();

			for (int i = 0; i < _master->getNumVars(); i++) {
				if (cstat[i] == IloCplex::BasisStatus::AtUpper) {
					rhs -= row[i];
					row[i] = -row[i];
				}
			}
				

			if (CGEnv::isFractional(rhs)) {
				
				double lb = CGEnv::FractionalPart(rhs);
				assert(lb >= 0.0);

				//CGConstraint* Cut = new CGConstraintGomory(_master->getEnv(), lb, IloInfinity);
				auto basisrow = _master->getInvBasisRow(k);
				CGConstraintGomory* Cut = new CGConstraintGomory(_master->getEnv(), lb, IloInfinity, CSTS, basisrow);
				bool isNonEmpty = false;
				IloNumArray coeffs(_master->getEnv().getImpl());
				IloNumVarArray vars(_master->getEnv().getImpl());


				for (int i = 0; i < _master->getNumVars(); i++) {
					CGNumVar* x = _master->getVar(i);

					//if (cstat[i] == IloCplex::BasisStatus::Basic) continue;

					double coef = CGEnv::FractionalPart(row[i]);
					assert(coef > -CGEnv::EPSILON);
					if (fabs(coef) > CGEnv::EPSILON)
						isNonEmpty = true;
					else coef = 0.0;
					if ( cstat[i] == IloCplex::BasisStatus::AtUpper) {
						lb -= coef;
						coef = -coef;
					}
					coeffs.add(coef);
					vars.add(*x);
				}
				Cut->setLinearCoefs(vars, coeffs);
				coeffs.end();
				vars.end();
				Cut->setLB(lb);

				if(!isNonEmpty)
					throw runtime_error("Empty Gomory Cut");

				//if (_master->acceptCut(*Cut, vals))
				//	CutsCollected.push_back({ Cut, 0.0 });
				//else
				//	delete Cut;

				double d = _master->getDistanceToCut(*Cut, vals);
				if (d > CGEnv::DELTA /*&& lb>CGEnv::EPSILON*/) {
					//Cuts.push_back(Cut);
					//cout << *Cut << " - " << d << endl;
					CutsCollected.push_back({ Cut, d });
				}
				else {
					delete Cut;
				}
			}

			cstat.end();
		}
	}

	//duals.end();

	delete[] bhead.first;
	delete[] bhead.second;

	//We sort and add the most violated Cuts
	sort(CutsCollected.begin(), CutsCollected.end(), [](pair<CGConstraint*, double>& a, pair<CGConstraint*, double>& b) {return a.second > b.second; });

	int maxcuts = 10;
	for (int i = 0; i < min(maxcuts, (int)CutsCollected.size()); i++) {
		Cuts.push_back(CutsCollected[i].first);
	}

	//We delete the remaining cuts
	for (int i = maxcuts; i < CutsCollected.size(); i++) {
		delete CutsCollected[i].first;
	}

	cout << "#Gomory Found: " << Cuts.size() - initSize << endl;

	return initSize < Cuts.size();
}

bool BBCuttingPlaneSolver::getMatchingCut(IloNumArray& vals , vector<CGConstraint*>& Cuts)
{

	if (CGEnv::isInteger(vals)) return false;


	//We get the unique values in vals
	vector<double> X;

	for (int k = 0; k < vals.getSize(); k++) /*if(!CGEnv::isInteger(vals[k]))*/{
		if (vals[k] > CGEnv::EPSILON && _master->isEdgeVar(k)) X.push_back(vals[k]);
	}
	sort(X.begin(), X.end(), [](double x1, double x2) {return x1 < x2; });

	//auto it = unique(X.begin(), X.end());
	//X.resize(distance(X.begin(), it));

	

    X.erase(std::unique(X.begin(), X.end(), [](double a, double b) { return std::abs(a - b) < 0.1; }), X.end());

	sort(X.begin(), X.end(), [](double x1, double x2) {return x1 > x2; });
	//X.clear();
	//X.push_back(1.0 - CGEnv::EPSILON);

	//We get the edges
	vector<CGEdge> E;
	for (int k = 0; k < vals.getSize(); k++) {
		if (vals[k] > CGEnv::EPSILON && _master->isEdgeVar(k)) {
			CGNumVar* x = _master->getVar(k);
			E.push_back({ x->geti(), x->getj(), vals[k] });
		}
	}
	sort(E.begin(), E.end(), [](CGEdge& a, CGEdge& b) {return a.w < b.w; });

	size_t initSize = Cuts.size();

	bool cutfound = false;
	for (auto threshold : X) if (!cutfound ) {
		vector<int> parent(_I->nnodes(), -1);
		vector<int> rank(_I->nnodes(), 1);

		vector<CGEdge> Et;
		copy_if(E.begin(), E.end(), back_inserter(Et), [threshold](CGEdge& e) {return e.w < threshold; });
		sort(Et.begin(), Et.end(), [](CGEdge& a, CGEdge& b) {return a.w < b.w; });

		vector<CGEdge> Tree;
		for (auto& currentedge : Et) {
			if (findparent(currentedge.i, parent) != findparent(currentedge.j, parent)) {
				unite(currentedge.i, currentedge.j, parent, rank);
				Tree.push_back(currentedge);

				//We get the comb
				vector<int> H;
				stack<int> Stack;
				Stack.push(currentedge.i);
				vector<bool> visited(_I->nnodes(), false);
				
				while (!Stack.empty()) {
					int i = Stack.top();
					Stack.pop();
					if (visited[i])
						continue;
					H.push_back(i);
					visited[i] = true;
					for (auto& e : Tree)
						if (e.i == i && !visited[e.j]) Stack.push(e.j);
						else if (e.j == i && !visited[e.i]) Stack.push(e.i);

				}

				//We get the arcs leaving H in the current Tree
				vector<bool> isCandidate(_I->nnodes(), true);
				for (auto i : H) isCandidate[i] = false;
				vector< CGEdge > DeltaH;
				copy_if(E.begin(), E.end(), back_inserter(DeltaH), [isCandidate,&threshold](CGEdge& e) {return isCandidate[e.i] != isCandidate[e.j] && e.w>=threshold; });
				//for (auto& e : E) //We consider all possible edges not only Et
				//	if ((isCandidate[e.i] && !isCandidate[e.j]) || (isCandidate[e.j] && !isCandidate[e.i]))
				//		DeltaH.push_back(e);

				if (DeltaH.empty()) continue;

				//We repair the teeth
				unordered_set< int > VT;
				for (auto& e : DeltaH) {
					VT.insert(e.i);
					VT.insert(e.j);
				}
				for (auto i : VT) {
					auto n = count_if(DeltaH.begin(), DeltaH.end(), [i](CGEdge& a) {return a.i == i || a.j == i; });
					//assert(n >= 1);
					if (n >= 2) {
						auto it = find(H.begin(), H.end(), i);
						if (it != H.end()) H.erase(it);
						else H.push_back(i);
						auto last = remove_if(DeltaH.begin(), DeltaH.end(), [i](CGEdge& a) {return a.i == i || a.j == i; });
						assert(last != DeltaH.end());
						DeltaH.erase(last, DeltaH.end());
					}
				}
				if (DeltaH.empty()) continue;


				//We get the Teeth
				sort(DeltaH.begin(), DeltaH.end(), [](CGEdge& a, CGEdge& b) {return a.w < b.w; });
				//vector< CGEdge > T;
				//T.push_back(DeltaH.back());
				//DeltaH.pop_back();
				//while (DeltaH.size() >= 2) {
				//	CGEdge e1 = DeltaH.back(); DeltaH.pop_back();
				//	CGEdge e2 = DeltaH.back(); DeltaH.pop_back();
				//	if (e1.w + e2.w - 1 <= CGEnv::EPSILON) break;
				//	else {
				//		T.push_back(e1);
				//		T.push_back(e2);
				//	}
				//}
				vector< CGEdge > T(DeltaH.begin(), DeltaH.end());

				////We repair the teeth
				////auto copyT(T);
				//unordered_set< int > VT;
				//for (auto& e : T) {
				//	VT.insert(e.i);
				//	VT.insert(e.j);
				//}
				//for (auto i : VT) {
				//	auto n = count_if(T.begin(), T.end(), [i](CGEdge& a) {return a.i == i || a.j == i; });
				//	//assert(n >= 1);
				//	if (n >= 2) {
				//		auto it = find(H.begin(), H.end(), i);
				//		if (it != H.end()) H.erase(it);
				//		else H.push_back(i);
				//		auto last = remove_if(T.begin(), T.end(), [i](CGEdge& a) {return a.i == i || a.j == i; });
				//		assert(last != T.end());
				//		T.erase(last, T.end());
				//	}
				//}

				auto c = count_if(H.begin(), H.end(), [currentedge](int i) {return currentedge.i == i || currentedge.j == i; });
				if (c ==0) continue;

				if (T.size() < 1) continue;
				if (T.size() % 2 == 0) {
					auto it = min_element(T.begin(), T.end(), [](CGEdge& a, CGEdge& b) {return a.w < b.w; });
					//We had the missig head to H
					CGEdge e = *it;
					if (find(H.begin(), H.end(), e.i) == H.end()) H.push_back(e.i);
					if (find(H.begin(), H.end(), e.j) == H.end()) H.push_back(e.j);
					T.erase(it);
				}
				if (T.size() >= 3) {
					
					CGConstraint* Cut = new CGConstraintMatching(_master->getEnv(), H, T);
					_master->setupConstraint(Cut);

					if (_master->acceptCut(*Cut, vals)) {
						Cuts.push_back(Cut);
						cutfound = true;
					}
					else
						delete Cut;
				}
				else {
					//We add the other head of edges in T to H
					//for (auto& e : T) {
					//	if (find(H.begin(), H.end(), e.i) == H.end()) H.push_back(e.i);
					//	if (find(H.begin(), H.end(), e.j) == H.end()) H.push_back(e.j);
					//}

					CGConstraint* cut = new CGConstraintConnectivity(_master->getEnv(), H);
					_master->setupConstraint(cut);
					if (_master->getDistanceToCut(*cut, vals) > CGEnv::EPSILON)
						Cuts.push_back(cut);
					else delete cut;
				}
				

			}
		}
	}

	cout << "#Matching Found: " << Cuts.size() - initSize << endl;

	return initSize < Cuts.size();
}

bool BBCuttingPlaneSolver::getChavtalCuts(IloNumArray& vals, vector<CGConstraint*>& Cuts)
{
	//Similar to Gomory but we check combination of pait of constraints divided by 2
	//{0, 1/2}-Chvatal-Gomory Cut
	Instance* I = _master->getEnv().getInstance();
	int n = I->nnodes();
	size_t initSize = Cuts.size();
	auto duals = _master->getDuals();
	for (int ki = 0; ki < _master->getNumConstraints(); ++ki) {
		for (int kj = ki + 1; kj < _master->getNumConstraints(); ++kj) {
			CGConstraint* Ri = _master->getConstraint(ki);
			CGConstraint* Rj = _master->getConstraint(kj);

			vector<CGConstraint*> R = { Ri, Rj };
			vector<double> coeffs = { 0.5, 0.5 };
			if (duals[ki] > CGEnv::EPSILON && duals[kj] > CGEnv::EPSILON)
				coeffs = { duals[ki], duals[kj] };

			double rhs = 0.0, lhs = 0.0;
			for (int c = 0; c < coeffs.size(); c++) {
				rhs += coeffs[c] * R[c]->getUB();
				lhs += coeffs[c] * R[c]->getLB();
			}
			assert(fabs(rhs - lhs) <= CGEnv::EPSILON);

			if (CGEnv::isInteger(lhs)) continue;

			lhs = floor(lhs);

			CGConstraint* Cut = new CGConstraintChvatal(_master->getEnv(), lhs, R, coeffs);

			bool isEmpty = true;
			for (int k = 0; k < vals.getSize(); ++k) {
				CGNumVar* x = _master->getVar(k);
				double coef = Cut->getCoefficient(x->geti(), x->getj());
				Cut->setLinearCoef(*x, coef);
				if (fabs(coef) > CGEnv::EPSILON) isEmpty = false;
			}

			if (!isEmpty && _master->acceptCut(*Cut, vals) > CGEnv::EPSILON) {
				double d = _master->acceptCut(*Cut, vals);

				Cuts.push_back(Cut);
				//cout << "Chvatal Cut: " << Cut->getLB() << " - " << Cut->getUB() << endl;
				break;
			}
			else
				delete Cut;
		}
	}
	duals.end();
	cout << "#Chvatal Found: " << Cuts.size() - initSize << endl;

	return initSize < Cuts.size();

}

bool BBCuttingPlaneSolver::getInfeasiblePathCuts(IloNumArray& vals, vector<CGConstraint*>& Cuts, double IncObj)
{
	//We get the edges in the solution
	vector<CGEdge> E;
	for (int k = 0; k < vals.getSize(); k++) {
		if (vals[k] > CGEnv::EPSILON && _master->isEdgeVar(k)) {
			CGNumVar* x = _master->getVar(k);
			E.push_back({ x->geti(), x->getj(), vals[k] });
		}
	}

	//We sort the edges by weight
	sort(E.begin(), E.end(), [&](CGEdge& a, CGEdge& b) {return a.w / _master->getLength(a.i, a.j) > b.w / _master->getLength(b.i, b.j); });

	//We greedily select an edge until the L of col obj is greater than the incumbent
	double L = 0.0;
	double sum = 0.0;
	vector<CGEdge> P;
	for (auto& e : E) {
		L += _master->getLength(e.i, e.j);
		sum += e.w;
		P.push_back(e);
		if (L >= IncObj) break;
	}

	if(L >= IncObj && sum > P.size()-1)
		cout << "Infeasible Path Cut: " << L << " - " << IncObj << " - " << sum << " - " << P.size() << endl;


	return false;
}

BCPSolver::BCPSolver(Instance* I) : Solver(I, "BCP")
{
	_master = new CGMasterProblem(I);
	_cuttingplanesolver = new BBCuttingPlaneSolver(_master, [this](IloNumArray& vals) { return applyHeuristic(vals); });
	_bestobjval = numeric_limits<double>::max();

}

BCPSolver::~BCPSolver()
{
	delete _cuttingplanesolver;
	delete _master;
	for (auto n : _BBNodes) delete n;
}

Solution* BCPSolver::recoversolution()
{
	if (_besttour.empty()) return nullptr;

	//We recover the best tour
	Solution* S = new Solution(_I);
	for (int i = 0; i < _besttour.size(); i++) {
		S->push_back(_besttour[i]);
	}
	return S;
}

void BCPSolver::addNode(BBNode* node)
{
	// We set the pending status
	node->setStatus(CGEnv::PendingNode);

	// We perform an insertion in decreasing order in _BBNodes
	auto it = _BBNodes.begin();
	while (it != _BBNodes.end() && (*it)->getObjVal() >= node->getObjVal()) {
		++it;
	}
	_BBNodes.insert(it, node);
	
}

bool BCPSolver::applyHeuristic(IloNumArray& vals)
{

	vector<int> T = extractTour(vals);
	return addFeasibleTour(T);

	//vector<int> tour;
	//vector<bool> visitednode(_I->nnodes(), false);
	//int n = _I->nnodes();


	////We get the first edge
	//CGEdge ebest = { -1, -1, -1 };
	//double valbest = -1;
	//for (int k = 0; k < vals.getSize(); k++) if(_master->getVar(k)->getType()==CGEnv::EdgeVar){
	//	int i = _master->getVar(k)->geti();
	//	int j = _master->getVar(k)->getj();
	//	double w = _I->travellingtime(_I->Node(i), _I->Node(j));
	//	if (vals[k] > valbest /*|| (fabs(vals[k] - valbest) < CGEnv::EPSILON && w < ebest.w)*/) {
	//		valbest = vals[k];
	//		ebest = { _master->getVar(k)->geti(), _master->getVar(k)->getj(), w };
	//	}
	//	else if (fabs(vals[k] - valbest) < CGEnv::EPSILON && RAND01()<0.5) {
	//		valbest = vals[k];
	//		ebest = { _master->getVar(k)->geti(), _master->getVar(k)->getj(), w };
	//	}
	//}
	//tour.push_back(ebest.i);
	//tour.push_back(ebest.j);
	//n -= 2;

	//while (n != 0)
	//{
	//	//we get the extremities of the tour
	//	int i = tour.front();
	//	int j = tour.back();
	//	visitednode[i] = true;
	//	visitednode[j] = true;

	//	//We get the best edge connected to one of the extremities
	//	CGEdge ebest = { -1, -1, DBL_MAX };
	//	double valbest = -1;
	//	for (int k = 0; k < vals.getSize(); k++) if(valbest < 1.0 - CGEnv::EPSILON) {
	//		CGNumVar* x = _master->getVar(k);
	//		if (x->getType() == CGEnv::SlackVar) continue;
	//		int ki = x->geti();
	//		int kj = x->getj();
	//		if ((!visitednode[ki] || !visitednode[kj])) {
	//			if (ki == i || kj == i || ki == j || kj == j) {
	//				double w = _I->travellingtime(_I->Node(ki), _I->Node(kj));
	//				if (vals[k] > valbest || (vals[k] == valbest && w < ebest.w)) {
	//					ebest = { ki, kj, w };
	//					valbest = vals[k];
	//				}
	//			}
	//		}
	//	}

	//	if (ebest.i == -1 || ebest.j == -1) {
	//		//We may fall here if there are fractional subtours
	//		//We get the best edge connected to one of the extremities
	//		ebest = { -1, -1, DBL_MAX };
	//		for (int ki = 0; ki < _I->nnodes(); ki++) {
	//			for (int kj = 0; kj < ki; kj++){
	//				if ((!visitednode[ki] || !visitednode[kj])) {
	//					if (ki == i || kj == i || ki == j || kj == j) {
	//						double w = _I->travellingtime(_I->Node(ki), _I->Node(kj));
	//						if (w < ebest.w) {
	//							ebest = { ki, kj, w };
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}

	//	if (ebest.i == -1 || ebest.j == -1)
	//		throw std::runtime_error("No edge found (Heuristic)");


	//	//We add the edge to the tour
	//	if (ebest.i == i) tour.insert(tour.begin(), ebest.j);
	//	else if (ebest.j == i) tour.insert(tour.begin(), ebest.i);
	//	else if (ebest.i == j) tour.push_back(ebest.j);
	//	else tour.push_back(ebest.i);

	//	n--;
	//}

	////We add the feasible tour
	//return addFeasibleTour(tour);
}

vector<int> BCPSolver::extractTour(IloNumArray& vals)
{
	//We first build a random path similar to GRASP and then we do a best insertion to insert missing nodes
	vector<int> tour;

	//We get the edges of the current solution
	vector<CGEdge> E;
	for (int k = 0; k < vals.getSize(); k++) {
		if (vals[k] > CGEnv::EPSILON && _master->isEdgeVar(k)) {
			CGNumVar* x = _master->getVar(k);
			E.push_back({ x->geti(), x->getj(), vals[k] });
		}
	}

	//We randomly select connected edge to construct the tour using a roulette wheel selection
	vector<int> visited(_I->nnodes(), false);
	int current = rand() % _I->nnodes();
	tour.push_back(current);
	visited[current] = true;
	while (!E.empty()) {
		//We get the extremities
		int i = tour.back();
		int j = tour.front();
		//We get the best edge
		CGEdge ebest = { -1, -1, -1 };
		for (auto& e : E) {
			if (e.i == i || e.j == i || e.i == j || e.j == j) {
				if (e.w > ebest.w || (fabs(e.w - ebest.w) < CGEnv::EPSILON && _master->getLength(e.i, e.j) < _master->getLength(ebest.i, ebest.j))) {
					ebest = e;
				}
			}
		}
		if (ebest.i == -1 || ebest.j == -1) break;

		if (ebest.i == i || ebest.j == i)
			tour.push_back(ebest.i == i ? ebest.j : ebest.i);
		else
			tour.insert(tour.begin(), ebest.i == j ? ebest.j : ebest.i);

		visited[tour.back()] = true;
		visited[tour.front()] = true;

		auto it = remove_if(E.begin(), E.end(), [visited](CGEdge& e) {return visited[e.i] && visited[e.j]; });
		E.erase(it, E.end());
	}

	//We complete the tour with a best insertion
	vector<int> missing;
	for (int i = 0; i < _I->nnodes(); i++) if (!visited[i]) missing.push_back(i);
	//while (!missing.empty()) {
	//	double best = DBL_MAX;
	//	int besti = -1;
	//	int bestj = -1;
	//	for (int i = 0; i < tour.size(); i++) {
	//		int k = (i + 1) % tour.size();
	//		double lik = _master->getLength(tour[i], tour[k]);
	//		for (int j = 0; j < missing.size(); j++) {
	//			double w = _master->getLength(tour[i], missing[j]);
	//			w += _master->getLength(missing[j], tour[k]);
	//			w -= lik;
	//			if (w < best) {
	//				best = w;
	//				besti = i;
	//				bestj = j;
	//			}
	//		}
	//	}
	//	tour.insert(tour.begin() + besti + 1, missing[bestj]);
	//	missing.erase(missing.begin() + bestj);
	//}

	//We do a best insertion of the missing node in a random order
	while (!missing.empty()) {
		//We get a random node
		int k = rand() % missing.size();
		int j = missing[k];
		missing.erase(missing.begin() + k);

		int besti = -1;
		double best = DBL_MAX;
		for (int i = 0; i < tour.size(); i++) {
			int l = (i + 1) % tour.size();
			double w = _master->getLength(tour[i], j) + _master->getLength(j, tour[l]) - _master->getLength(tour[i], tour[l]);
			if (w < best) {
				best = w;
				besti = i;
			}
		}

		tour.insert(tour.begin() + besti + 1, j);
	}




	return tour;
}

bool BCPSolver::addFeasibleTour(vector<int>& tour)
{
	double length = _I->length(tour);
	if (length < _bestobjval) {
		_bestobjval = length;
		_besttour = tour;
		cout << "New Incumbent: " << _bestobjval << endl;

		//We remove bounded nodes
		if (!_BBNodes.empty()) 
			cout << "Obj Range: [" << _BBNodes.front()->getObjVal() << ", " << _BBNodes.back()->getObjVal() <<"]" << endl;
		auto it = _BBNodes.begin();
		while (it != _BBNodes.end()) {
			if ((*it)->getObjVal() > _bestobjval) {
				delete* it;
				it = _BBNodes.erase(it);
			}
			else break;
		}

		return true;
	}
	else return false;
}

vector<CGNumVar*> BCPSolver::getBranchingCandidates(IloNumArray& vals)
{
	//Breaking up suptour of three or list of fractional variables
	vector<CGNumVar*> candidates;

	vector<CGNumVar*> fractionalvars;
	vector< vector<CGNumVar*> > connected(_I->nnodes());
	for (int k = 0; k < vals.getSize(); k++) {
		if (_master->isEdgeVar(k) && vals[k]>CGEnv::EPSILON/*CGEnv::isFractional(vals[k])*/) {
			if(CGEnv::isFractional(vals[k])) fractionalvars.push_back(_master->getVar(k));
			connected[_master->getVar(k)->geti()].push_back(_master->getVar(k));
			connected[_master->getVar(k)->getj()].push_back(_master->getVar(k));
		}
	}

	return fractionalvars;

	//Now we identify all subtours of three
	for (int i = 0; i < _I->nnodes(); i++) {
		for (auto x : connected[i]) for (auto y : connected[i]) if (x != y) {
			//we get the extremities other than i
			int j = x->geti() == i ? x->getj() : x->geti();
			int k = y->geti() == i ? y->getj() : y->geti();
			//We then check in connected[j] if there is an edge between j and k
			for (auto z : connected[j]) {
				if ((z->geti() == k || z->getj() == k)) {
					if (CGEnv::isFractional(vals[x->getId()])) candidates.push_back(x);
					if (CGEnv::isFractional(vals[y->getId()])) candidates.push_back(y);
					if (CGEnv::isFractional(vals[z->getId()])) candidates.push_back(z);
				}
			}
		}
	}

	//return fractionalvars;

	if (candidates.empty()) return fractionalvars;
	else return candidates;
	

}

double BBFlowGraph::getFlow(const vector<int>& shore1, const vector<int>& shore2)
{
	//We get the total flow between the shores
	double flow = 0.0;
	for (int i : shore1) {
		for (int k : _adj[i]) {
			FlowEdge* e = _edges[k];
			if (shore2.end() != find(shore2.begin(), shore2.end(), e->j)) {
				flow += e->w;
			}
		}
	}
	return flow;
}

double BBFlowGraph::getLeavingFlow(const vector<int>& shore)
{
	//We get the total flow leaving the shore
	double flow = 0.0;
	for (int i : shore) {
		for (int k : _adj[i]) {
			FlowEdge* e = _edges[k];
			if (shore.end() == find(shore.begin(), shore.end(), e->j)) {
				flow += e->w;
			}
		}
	}
	return flow;
}

double BBFlowGraph::getFlow(vector<FlowEdge*>& P)
{
	//We get the flow
	double flow = DBL_MAX;
	for (auto e : P) {
		double f = residualflow(e);
		if (f < flow) flow = f;
	}
	return flow;
}

vector<FlowEdge*> BBFlowGraph::getPath()
{
	vector<FlowEdge*> Path;

	fill(_pred.begin(), _pred.end(), -2);

	for (auto t : _sinks) _pred[t] = -3;

	vector<int> Q;
	for (auto s : _sources) {
		Q.push_back(s);
		_pred[s] = -1;
	}


	while (!Q.empty()) {
		int i = Q.back();
		Q.pop_back();

		for (int k : _adj[i]){
			FlowEdge* e = _edges[k];
			int j = e->j;

			if (residualflow(e) < CGEnv::EPSILON) continue;

			//We check if we reach one of the targets
			if (_pred[j] == -3) {
				_pred[j] = e->id;
				Q.clear();
				for (int k = _pred[j]; k >= 0; k = _pred[_edges[k]->i]) Path.push_back(_edges[k]);
				reverse(Path.begin(), Path.end());
				break;
			}
			else if (_pred[j] == -2 /*&& e->w - e->f > CGEnv::EPSILON*/ ) {
				_pred[j] = e->id;
				Q.push_back(j);
			}
		}
	}


	return Path;

}

double BBFlowGraph::maxFlow()
{
	resetFlow();
	double flow = 0.0;
	vector<FlowEdge*> P = getPath();
	while (!P.empty() && flow < 2) {
		double f = getFlow(P);
		flow += f;
		for (auto e : P) {
			e->f -= f;
			e->r->f += f;
		}
		P = getPath();
	}

	return flow;

}

double BBFlowGraph::getShores(vector<int>& shore1, vector<int>& shore2, bool strict)
{

	//We get the max flow
	double flow = maxFlow();
	//strict = false;
	//We get the shores
	if(!strict){
		for (int i = 0; i < _nnodes; i++) {
			if (_pred[i] >= -1) shore1.push_back(i);
			else shore2.push_back(i);
		}
	}
	else{
		vector<int> rank(_nnodes, max(2, _nnodes));
		for (int s : _sources) rank[s] = -1;
		for (auto t : _sinks) rank[t] = -2;

		//We propagate the ranks
		auto it = find_if(rank.begin(), rank.end(), [&](int r) { return r < 0; });
		while (it != rank.end()) {
			int i = std::distance(rank.begin(), it);
			assert(rank[i] < 0);
			int ri = -rank[i];
			rank[i] = ri;
			for (int k : _adj[i]) {
				FlowEdge* e = _edges[k];
				if (residualflow(e) < CGEnv::EPSILON) continue;
				int j = (e->i == i) ? e->j : e->i;
				int rj = fabs(rank[j]);
				if (rj > 0 && rj > ri) {
					rank[j] = -ri;
				}
			}
			it = find_if(rank.begin(), rank.end(), [&](int r) { return r < 0; });
		}

		for (int i = 0; i < _nnodes; i++) {
			if (rank[i] == 1) shore1.push_back(i);
			else if (rank[i] == 2) shore2.push_back(i);
		}
	}

	//vector<int> s1, s2;
	//for (int i = 0; i < _nnodes; i++) {
	//	if (_pred[i] >= -1) s1.push_back(i);
	//	else s2.push_back(i);
	//}

#if DEBUG
	//We recompute the flow to verify the correctness
	double flow1 = getLeavingFlow(shore1);
	double flow2 = getLeavingFlow(shore2);
	assert(fabs(flow - flow1) < CGEnv::EPSILON);
	assert(fabs(flow - flow2) < CGEnv::EPSILON);
#endif

	return flow;
}

void BBFlowGraph::exportGraph(string name)
{
	//We export the graph in a text file for graphviz visualization
	ofstream file(name);
	file << "digraph G {" << endl;
	for (int i = 0; i < _nnodes; i++) {
		file << i << " [label=\"" << i << "\"];" << endl;
	}
	for (int i = 0; i < _nnodes; i++) {
		for (int k : _adj[i]) {
			FlowEdge* e = _edges[k];
			file << e->i << " -> " << e->j << " [label=\"" << e->w << "/" << e->f << "\"];" << endl;
		}
	}
	file << "}" << endl;
	file.close();
}

bool BBCuttingPlaneSolver::getFlowCuts(IloNumArray& vals, vector<CGConstraint*>& Cuts)
{
	//We init the flow graph
	BBFlowGraph G(_I->nnodes());
	for (int k = 0; k < _master->getNumVars(); k++) if (_master->isEdgeVar(k)){
		double val = vals[k];
		if (val > CGEnv::EPSILON) {
			
			CGNumVar* x = _master->getVar(k);
			G.addEdge(x->geti(), x->getj(), val, val);
		}
		
	}
	//for (int k = 0; k < _master->getNumVars(); k++) if (!_master->isEdgeVar(k))
	//	cout << vals[k] << " ";
	//cout << endl;

	//We get the potential sources
	vector<int> Sources;
	int s = rand() % _I->nnodes();
	Sources.push_back(s);
	//for (int i = 0; i < _I->nnodes(); i++) Sources.push_back(i);

	size_t initSize = Cuts.size();

	while (!Sources.empty()) {
		
		int s = Sources.back();
		Sources.pop_back();

		G.resetSources();
		G.addSource(s);

		//We create the potential sinks
		vector<int> Sinks;
		for (int i = 0; i < _I->nnodes(); i++) if (i != s) Sinks.push_back(i);

		while (!Sinks.empty()) {
			int t = Sinks.back();
			Sinks.pop_back();

			G.resetSinks();
			G.addSink(t);

			//We get the shores
			vector<int> shore1, shore2;	
			double flow = G.getShores(shore1, shore2, true);

			
			

			size_t ncuts = Cuts.size();
			if (flow < 2.0 - CGEnv::EPSILON) {
				bool cutfound = false;
				//We add the connectivity cuts
				for (auto& shore : { shore1, shore2 }) {
					CGConstraint* cut = new CGConstraintConnectivity(_master->getEnv(), shore);
					_master->setupConstraint(cut);
					if (_master->acceptCut(*cut, vals)) {
						Cuts.push_back(cut);
						cutfound = true;
					}
					else delete cut;
				}

				if (!cutfound) {
					//We add the subtour elimination cuts
					for (auto& shore : { shore1, shore2 }) if (shore.size() > 1) {
						CGConstraint* cut = new CGConstraintSubtourElimination(_master->getEnv(), shore);
						_master->setupConstraint(cut);
						//cout << "Flow Cut: " << "OutFlow:" << " = " << G.getLeavingFlow(shore) << endl;
						//cout << "Flow Cut: " << "InFlow:" << " = " << G.getFlow(shore, shore) << endl;
						//cout << _master->getDistanceToCut(*cut, vals) << endl;
						if (_master->acceptCut(*cut, vals)) {
							Cuts.push_back(cut);
						}
						else delete cut;
					}
				}
				
				
			}

			if (ncuts<Cuts.size()) {
				vector<bool> inshore(_I->nnodes(), false);
				for (auto& shore : { shore1, shore2 })
					for (int i : shore) inshore[i] = true;



				//We remove the nodes in the shore from the sinks
				auto itsinks = remove_if(Sinks.begin(), Sinks.end(), [&](int i) { return inshore[i]; });
				Sinks.erase(itsinks, Sinks.end());
			}
			
				
		}
	}

	//if (Cuts.size() - initSize ==0 ) {
	//	G.exportGraph("flowgraph.dot");
	//	exit(222);
	//}

	//assert(!CGEnv::isInteger(vals) || Cuts.size() - initSize > 0);

	cout << "#Flow Found: " << Cuts.size() - initSize << endl; 

	return initSize < Cuts.size();
}

void BBCuttingPlaneSolver::solve(BBNode* N)
{
	auto start = chrono::steady_clock::now();
	bool hasNewCuts = true;
	int iter = 0;
	double bestobj = 0.0;
	while (iter++<10) {

		//We solve the master problem
		CGEnv::SolverStatus status = _master->solve(N);
		

		if (status == CGEnv::Infeasible) {
			throw runtime_error("Infeasible problem at the root");
		}

		double objval = _master->getObjVal();

		double gap = (objval - bestobj) / max(1.0, fabs(objval));
		if (gap >= 1e-4 ) {
			iter = 0;
			bestobj = objval;
		}
		

		//We get the solution
		IloNumArray vals = _master->getVarVals();

		//if (heuristicFunction != nullptr) heuristicFunction(vals);

		//if (CGEnv::isInteger(vals)) break;

		vector<CGConstraint*> Cuts;

		hasNewCuts = false;

		//We get the Chvatal cuts
		//hasNewCuts = hasNewCuts || getChavtalCuts(vals, Cuts);
		
		
		
		//We get the matching cuts
		//hasNewCuts = hasNewCuts || getMatchingCut(vals, Cuts) || hasNewCuts;
		//hasNewCuts = Cuts.size() >= 10;

		//We get the flow cuts
		hasNewCuts = /*hasNewCuts ||*/ getFlowCuts(vals, Cuts) || hasNewCuts;

		//if (iter >= 5 ) hasNewCuts = getFlowCuts(vals, Cuts) || hasNewCuts;

		
		//hasNewCuts = Cuts.size() >= 10;
		//hasNewCuts = /*hasNewCuts ||*/ getGomoryCuts(vals, Cuts) || hasNewCuts;

		//We get the gomory cuts: very slow
		//if(iter>=5) hasNewCuts = getGomoryCuts(vals, Cuts) || hasNewCuts;
		hasNewCuts = hasNewCuts || getGomoryCuts(vals, Cuts) || hasNewCuts;

		

		//We add the cut to the model
		_master->addConstraints(Cuts);
		

		//WE get the number of added constraints to print the log
		int ncuts = Cuts.size();

		cout << "CP: Obj=" << objval << ", #addedcsts = " << ncuts << ", #csts = " << _master->getNumConstraints() << endl;
	
		vals.end();
	}


	auto end = chrono::steady_clock::now();

	cout << "Elapsed time: " << chrono::duration_cast<chrono::seconds>(end - start).count() << " s" << endl;


}
