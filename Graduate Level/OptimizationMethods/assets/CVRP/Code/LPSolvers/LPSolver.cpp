#include "LPSolver.hpp"

NumVar* LPSolver::getVar(FlowEdge* e)
{
	int i = e->geti();
	int j = e->getj();
	auto Vars = _env.getVars();

	for (NumVar* X : Vars) {
		if (X->getType() == Solver::EdgeVar && X->geti() == i && X->getj() == j) return X;
	}
	return nullptr;
}

void LPSolver::convertToIP() {
	if (!_isIP) _env.convertToIP();
	_isIP = true;
}

void LPSolver::convertToLP() {
	if (_isIP) _env.convertToLP();
	_isIP = false;
}

LPCplexSolver::LPCplexSolver(Instance* I) : LPSolver(I)  {


	int origin = _I->InitialDepotId();
	int destination = _I->FinalDepotId();

	//The variables xij
	for (int i = 0; i < _I->nnodes(); i++) {
		for (int j = 0; j < _I->nnodes(); j++) {
			if (i != j && i!=destination && j!=origin ) {
				NumVar* X = new NumVar(_env, i, j);
				_env.addVar(X);
			}
		}
	}

	//The variable q
	for (int i = 0; i < _I->nnodes(); ++i) {
		node* N = _I->Node(i);
		NumVar* q = new NumVar(_env, i);
		q->setLB(N->demand());
		q->setUB(_I->VehicleCap());
		if (N->isInitialDepot()) q->setUB(0);
		if (N->isFinalDepot()) 
			q->setLB(_I->VehicleCap());
		_env.addVar(q);
	}

	

	//Objective function
	_env.setObjective([&](NumVar* X) -> double {
		if (X->getType() == Solver::EdgeVar) {
			int i = X->geti();
			int j = X->getj();
			return _I->cost(i, j);
		}
		else return 0.0;
		});

	//for (int k = 0; k < _env.nvars(); k++) {
	//	NumVar* X = _env.getVar(k);
	//	if (X->getType() == Solver::EdgeVar) {
	//		int i = X->geti();
	//		int j = X->getj();
	//		_env.setObjVarCoef(X, _I->cost(i, j));
	//	}
	//}
	

	//C1: Depots constraints	
	Constraint* outdepot = new ConstraintNodeDegreeOut(_env, origin, _I->FleetSize());
	Constraint* indepot = new ConstraintNodeDegreeIn(_env, destination, _I->FleetSize());

	_env.addConstraint(outdepot);
	_env.addConstraint(indepot);


	//C2: Node degree constraints
	for (int i = 0; i < _I->nnodes(); i++) if(i!=origin && i!=destination){
		Constraint* degout = new ConstraintNodeDegreeOut(_env, i);
		Constraint* degin = new ConstraintNodeDegreeIn(_env, i);

		_env.addConstraint(degout);
		_env.addConstraint(degin);
	}

	//C3: Load propagation constraints
	auto Vars = _env.getVars();
	for (NumVar* X : Vars) {
		if (X->getType() == Solver::EdgeVar) {
			int i = X->geti();
			int j = X->getj();

			Constraint* loadprop = new ConstraintLoadPropagation(_env, i, j);
			_env.addConstraint(loadprop);

		}
	}

	_env.setTimLim(_timlim);
	_env.setThreads(_threads);
	_env.setMinGap(_mingap);
    // Set CPLEX to use dual simplex
    //_env.getCplex().setParam(IloCplex::RootAlg, IloCplex::Dual);
	//_env.exportModel("LP.lp");
	enableOutput(false);
}




void CuttingPlaneSolver::solvemethod(Solution* S) {
	_env.setTimLim(_timlim);
	_env.setThreads(_threads);
	_env.setMinGap(_mingap);

	_env.DebugCheckConstraints(_env.getConstraints());


	// 1. Solve initial LP relaxation
	_status = _env.runsolver();
	if (_isIP) return;
	if (_runonce && _hasrun) return;
	_hasrun = true;
	
	// 2. Cutting plane loop
	bool cutsAdded = true;
	int round = 0;
	const int maxNoImprove = 10;
	int noImproveCount = 0;
	double lastObjVal = _env.getObjVal();

	std::cout << "[CuttingPlane] Initial LP bound: " << lastObjVal << std::endl;

	while (cutsAdded && _status == Solver::SolverStatus::Optimal && noImproveCount < maxNoImprove && timeLeft()>EPSILON) {
		cutsAdded = false;

		std::vector<Constraint*> Cuts;

		// 2a. Check for violated cuts
		std::vector<Constraint*> SECCuts = separateSECCuts();
		Cuts.insert(Cuts.end(), SECCuts.begin(), SECCuts.end());

		if(Cuts.empty()) {
			std::vector<Constraint*> GomoryCuts = separateGomoryCuts(0.1 * _env.getNumConstraints()); // Adjust the size as needed
			Cuts.insert(Cuts.end(), GomoryCuts.begin(), GomoryCuts.end());
		}

		int cutsThisRound = Cuts.size();
		cutsAdded = !Cuts.empty();

		// 2b. Add cuts if found
		_env.addConstraints(Cuts);
		
		// 2c. Re-solve if cuts were added
		if (cutsAdded) {
			_status = _env.runsolver();
			double currObjVal = _env.getObjVal();
			double percentImprovement = 0.0;
			if (std::abs(lastObjVal) > EPSILON) {
				percentImprovement = std::abs(lastObjVal - currObjVal) / std::max(std::abs(lastObjVal), std::abs(currObjVal));
			}
			// Log the number of cuts added and the current bound
			std::cout << "[CuttingPlane] Round " << ++round
				<< ": Cuts added = " << cutsThisRound
				<< ", Current LP bound = " << currObjVal 
				<< ", Improvement = " << percentImprovement
				<< ", Time = " << _timlim - timeLeft() << " s."
				<< std::endl;

            
            if (percentImprovement < 0.01) {
				noImproveCount++;
			} else {
				noImproveCount = 0;
			}
			lastObjVal = currObjVal;

			if (_status != Solver::SolverStatus::Optimal) {
				_env.exportModel("infeasible.lp");
				throw std::runtime_error("[CuttingPlane] LP became infeasible or unbounded after adding cuts.");
			}
		}
		else {
			std::cout << "[CuttingPlane] No cuts found in round " << ++round << ". Terminating." << std::endl;
		}
	}

	_env.DebugCheckConstraints(_env.getConstraints());

	Solution* Sol = recoversolution();
	if(Sol->isInteger()) cout << "[CuttingPlane] Found integer solution with objective value: " << Sol->fitness() << endl;
	else cout << "[CuttingPlane] Found fractional solution with objective value: " << Sol->fitness() << endl;
	delete Sol;


}	

std::vector<NumVar*> CuttingPlaneSolver::getsortedVars()
{
	auto names = _env.getSortedColNames();

	auto Vars = _env.getVars();

	assert(Vars.size() == names.size());

	// Store variables in Vars into sortedVars according to their name order in names
	std::vector<NumVar*> sortedVars;

	std::map<string, int> mapnames;
	for(int k=0; k<Vars.size(); k++) {
		mapnames[Vars[k]->getName()] = k;
	}

	for (const auto& varName : names) {
		int k = mapnames[varName];
		if (k < 0 || k >= Vars.size()) {
			throw std::runtime_error("Variable index out of bounds: " + varName);
		}
		sortedVars.push_back(Vars[k]);
	}

	return sortedVars;
}

std::vector<Constraint*> CuttingPlaneSolver::separateGomoryCuts(size_t SizeMax) {
	std::vector<Constraint*> cuts;

	int nRows = _env.getNumConstraints();
	int nVars = _env.nvars();

	auto Vars = getsortedVars();

	assert(Vars.size() == nVars);


	IloCplex::BasisStatusArray basisstatus = _env.getBasisStatus(Vars);
	auto bhead = _env.getBHead();

	std::vector<Constraint*> Constraints = _env.getConstraints();

	auto xVals = _env.getVarVals(Vars);

	//Filtering the constraints to find Gomory cuts
	vector<int> cstIds; //Store the constraint IDs for Gomory cuts
	vector<double> scores(nRows, 0.0); //Scores for the cuts
	for (int k = 0; k < nRows; ++k) {
		if (bhead.first[k] < 0) continue; // Skip non-basic rows
		assert(bhead.first[k] < nVars);
		NumVar* BasicVar = Vars[bhead.first[k]];
		double BasicVal = bhead.second[k];
		assert(basisstatus[bhead.first[k]] == IloCplex::BasisStatus::Basic);
		assert(fabs(BasicVal - xVals[bhead.first[k]]) < EPSILON);
		if (ORUtils::isInteger(BasicVal) || BasicVar->getType() != VarType::EdgeVar) continue;
		cstIds.push_back(k);
		scores[k] = _I->travellingtime(BasicVar->geti(), BasicVar->getj())*BasicVal; //Score based on travelling time
	}

	sort(cstIds.begin(), cstIds.end(), [&scores](int a, int b) {
		return scores[a] < scores[b]; // Sort by scores
		});


	//Gomory cuts generation
	for (int k : cstIds) {

		NumVar* BasicVar = Vars[bhead.first[k]];
		double BasicVal = bhead.second[k];

		double rhs = _env.getInvRHS(k);

		vector<double> tableauRow = _env.getInvRow(k);

		int nzcount = 0;
		for (int i = 0; i < nVars; ++i)
			if (fabs(tableauRow[i]) > EPSILON) ++nzcount;

		if (nzcount > nVars * 0.5) continue; // Skip overly dense rows

		//Adjust the tableau row and rhs to account for variables at their bounds
		std::vector<double> row;
		std::vector<double> col;
		bool isFractional = false; //Flag to check if the cut is fractional
		for (int i = 0; i < nVars; ++i) {
			NumVar* x = Vars[i];
			double xVal = xVals[i];
			//if (basisstatus[i] == IloCplex::BasisStatus::AtUpper) {
			//	row.push_back(tableauRow[i]);
			//	col.push_back(xVal);
			//	tableauRow[i] = 0.0;
			//}
			////else if (basisstatus[i] == IloCplex::BasisStatus::AtLower && x->getLB() > EPSILON) {
			////	row.push_back(tableauRow[i]);
			////	col.push_back(xVal);
			////	tableauRow[i] = 0.0;
			////}
			//else if (basisstatus[i] == IloCplex::BasisStatus::Basic) {
			//	if (x == BasicVar) assert(fabs(tableauRow[i] - 1.0) < EPSILON); //Basic variable should be 1
			//	else assert(fabs(tableauRow[i]) < EPSILON); //Other basic variables should be zero
			//}

			if(!ORUtils::isInteger(xVal)) {
				isFractional = true; //If any variable is fractional, the cut is fractional
			}
		}
		rhs = rhs - ORUtils::dotproduct(row, col); //Adjust rhs based on the tableau row and variable values




		double f0 = ORUtils::FractionalPart(rhs);

		if (f0 > 0.1 && f0<0.9 && isFractional) { 
			//cout << "f0: " << f0 << endl;
			std::vector<double> InvBasisRow = _env.getInvBasisRow(k);

			ConstraintGomory* cut = new ConstraintGomory(_env, f0, IloInfinity, Constraints, InvBasisRow);

			//We manually set the coefficients of the variables at their bounds for efficiency
			bool isempty = true;
			IloNumArray tableauRowArray(_env);
			IloNumVarArray tableauRowVars(_env);
			for (int i = 0; i < nVars; ++i) {
				NumVar* x = Vars[i];

				double coef =  ORUtils::FractionalPart(tableauRow[i]);
				
				if (fabs(coef) < EPSILON || coef>=1.0-EPSILON) continue; //Skip if the coefficient is zero


				tableauRowArray.add(coef);
				tableauRowVars.add(*x);
				assert(coef > 0 && coef < 1.0 - EPSILON);
				if (coef > EPSILON) isempty = false;
			}
			if(!isempty) cut->setLinearCoefs(tableauRowVars, tableauRowArray);

			tableauRowArray.end();
			tableauRowVars.end();

			if (!isempty && _env.acceptCut(*cut)) {
				cuts.push_back(cut);
			}
			else 
				delete cut;

		}

		if (cuts.size() >= SizeMax) break; //Stop if we reached the maximum number of cuts
	}

	xVals.end();
	basisstatus.end();
	delete[] bhead.first;
	delete[] bhead.second;

	return cuts;
}

std::vector<Constraint*> CuttingPlaneSolver::separateSECCuts(size_t SizeMax)
{
	std::vector<Constraint*> Cuts;

	Solution* S = recoversolution();
	assert(S != NULL);
	std::vector<bool> visited(_I->nnodes(), false);
	for (int i = 0; i < _I->nnodes(); ++i) {
		if (i == _I->InitialDepotId() || i == _I->FinalDepotId() || visited[i]) continue; //Skip depots

		std::vector<int> source;
		std::vector<int> sink;
		Constraint* cut = NULL;
		double flow = 0.0;
		std::pair< std::vector<int>, std::vector<int> > Shores;

		//i to destination
		source.push_back(i);
		sink.push_back(_I->FinalDepotId());

		Shores = S->maxFlow(source, sink);
		flow = S->getFlow(Shores.first, Shores.second);

		double rid = std::ceil(_env.getInstance()->demand(Shores.first) / _env.getInstance()->VehicleCap());
		if (flow < rid - EPSILON) {
			cut = new ConstraintConnectivity(_env, Shores.first, ConstraintConnectivity::direction::Out);
			//cut = new ConstraintSubtourElimination(_env, Shores.first);
			if (_env.acceptCut(*cut)) {
				Cuts.push_back(cut);
				for (auto id : Shores.first) {
					visited[id] = true; //Mark the nodes in the cut as visited
				}
			}
			else {
				delete cut;
				cut = NULL;
			}
		}

		

		source.clear();
		sink.clear();
		//if (cut != NULL) break;

		if (cut != NULL) continue; //If we already found a cut, skip the rest

		//continue; //for debugging purposes

		//origin to i
		source.push_back(_I->InitialDepotId());
		sink.push_back(i);
		Shores = S->maxFlow(source, sink);
		flow = S->getFlow(Shores.first, Shores.second);

		double roi = std::ceil(_env.getInstance()->demand(Shores.second) / _env.getInstance()->VehicleCap());;
		if (flow < roi - EPSILON) {
			cut = new ConstraintConnectivity(_env, Shores.second, ConstraintConnectivity::direction::In);
			//cut = new ConstraintSubtourElimination(_env, Shores.second);
			if (_env.acceptCut(*cut)) {
				Cuts.push_back(cut);
				for (auto id : Shores.second) {
					visited[id] = true; //Mark the nodes in the cut as visited
				}
			}
			else {
				delete cut;
				cut = NULL;
			}
		}




	}

	//if(Cuts.empty()) {
	//	S->exportGraph("solution.dot");
	//}

	delete S;

	return Cuts;
}
