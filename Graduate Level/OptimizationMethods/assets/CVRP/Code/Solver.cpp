#include "Solver.hpp"

unsigned int NumVar::_index = 0;
const double Env::DELTA = 0.001;

unsigned int Constraint::_index = 0;

inline string Solver::PadNumber(int num, int w, char c)
{
    std::ostringstream ss;
    ss << std::setw(w) << std::setfill(c) << num;
    return ss.str();
}

inline string Solver::RedString(string s, int w) {
    if (s.size() > w) {
        string f = "...";
        int ww = (int)floor((w - f.size()) / 2.0);
        return s.substr(0, ww) + f + s.substr(s.size() - ww, s.size());
    }
    else return s;
}


Solver::Solver(Instance* I, string name)
{
	_I = I;
	_name = name;
	_start=_end= chrono::system_clock::now();
	_timlim = 60;
	_mingap = 0.00;
    _threads = 1;
}


void Solver::setparam(param_t PARAM, double value)
{
	if (PARAM == TimLim) _timlim = value;
	else if (PARAM == Gap) _mingap = value;
    else if (PARAM == Threads) _threads = (int)fmax(1.0,value);
}

void Solver::solve(Solution* S)
{
	startTimer();

	solvemethod(S);

	stopTimer();
}


void Solver::save(string outputfile)
{
    if (outputfile == "None")
        outputfile = "Output/" + _name + "_" + to_string((int)_timlim) + ".txt";

    ofstream output;

    int width = 20;
    output.open(outputfile.c_str(), ios::app);

    output.setf(ios::left);
    output.precision(2);
    output.setf(ios::fixed, ios::floatfield);
    
    output << " ";
    if (output.tellp() == 1) {
        output << setw(width) << "Date" << setw(width) << "Time" << setw(40) << "Instance" << setw(width) << "Solver" << setw(width) << "Fitness" << setw(width) << "Gap" << setw(width) << "CPU_time" << setw(width) << "Feas" ;
        output << setw(width) << "#_nodes"<< endl;
        output << " ";
    }

    Solution* S = recoversolution();
    string fitness = "n/a";
    string feas = "n/a";
    if (S != NULL) {
        fitness = to_string(S->fitness());
        if (S->isFeasible(false)) feas = "yes";
        else feas = "no";
    }

    std::time_t end_time = std::chrono::system_clock::to_time_t(_end);
    std::tm now;
    localtime_s(&now, &end_time);

    string date = PadNumber(now.tm_mday) + "/" + PadNumber(now.tm_mon + 1) + "/" + to_string(now.tm_year + 1900);
    string time = PadNumber(now.tm_hour) + ":" + PadNumber(now.tm_min);

    output << setw(20) << date << setw(20) << time << setw(40) << RedString(_I->name(),40) << setw(width) << RedString(_name,width) << setw(width) << fitness << setw(width) << gap() << setw(width) << cpuTime() << setw(width) << feas ;
    output << setw(width) << _I->nnodes();

    output << endl;
    
    delete S;

    output.close();
}

string NumVar::nameEdgeGen(int i, int j) {
    return "X" + to_string(i) + "_" + to_string(j);
}

void NumVar::consolidate() {
    if (_type == Solver::UndEdgeVar) {
        if (_i < _j) swap(_i, _j);
        setName(NumVar::nameEdgeGen(_i, _j).c_str());
    }
    else if (_type == Solver::EdgeVar) {
        setName(NumVar::nameEdgeGen(_i, _j).c_str());
    }
    else if(_type == Solver::SlackVar){
        setName(NumVar::nameSlackGen(_id).c_str());
        _i = _j = _id; //Unique pairs
    }
    else if (_type == Solver::CapVar) {
        string s = "q_" + to_string(_i);
        setName(s.c_str());
    }
    else {
        string s = "x_" + to_string(_id);
        setName(s.c_str());
    }
}


   
void Constraint::addVars(vector<NumVar*>& X) {
    IloNumArray coefs(*_env);
    IloNumVarArray vars(*_env);
    for (auto& x : X) {
        double coef = this->getCoefficient(x);
        if (fabs(coef) > EPSILON) {
            coefs.add(coef);
            vars.add(*x);
        }
        
    }
    this->setLinearCoefs(vars, coefs);
    vars.end();
    coefs.end();
}

double Constraint::retrieveCoefficient(NumVar* x)
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
pair<NumVar*, double> Constraint::toStandarForm()
{
    //Transform the constraint in standard form by adding the slacks variables
    double lb = this->getLB();
    double ub = this->getUB();
    if (fabs(ub - lb) > EPSILON) {
        if (lb <= -IloInfinity && ub >= IloInfinity)
            throw runtime_error("Invalid constraint bounds");
        if(lb>-IloInfinity && ub < IloInfinity)
            throw runtime_error("Constraint must be separated for standard form");
        if (lb > -IloInfinity) {
            _SlackVariable = new NumVar(*_env);
            _slackcoef = -1.0;
            this->setLinearCoef(*_SlackVariable, _slackcoef);
            _env->addVar(_SlackVariable);
            this->setLB(lb);
            this->setUB(lb);
            
        }
        else if (ub < IloInfinity) {
            _SlackVariable = new NumVar(*_env);
            _slackcoef = 1.0;
            this->setLinearCoef(*_SlackVariable, 1.0);
            _env->addVar(_SlackVariable);
            this->setLB(ub);
            this->setUB(ub);
        }
        
    }
    return make_pair(_SlackVariable, _slackcoef);
}

bool Constraint::isValidUnderConstraint(const std::vector<NumVar*>& Vars, const std::vector<double>& Vals)
{
	double eval = 0.0;
    for (size_t k = 0; k < Vars.size(); ++k) {
        NumVar* var = Vars[k];
		Solver::VarType type = var->getType();
		if (type != Solver::EdgeVar && type != Solver::CapVar && type != Solver::SlackVar) continue; // Ignore non-edge variables

		double coef = this->getCoefficient(var);
		eval += coef * Vals[k];   
	}

	if (_SlackVariable == NULL) 
	    return eval >= this->getLB() - EPSILON && eval <= this->getUB() + EPSILON;
	else if (_slackcoef>0.5) // Positive slack variable
        return eval <= this->getUB() + EPSILON;
    else // Negative slack variable
		return eval >= this->getLB() - EPSILON;
}

void Constraint::print(const std::vector<NumVar*>& Vars, const std::vector<double>& Vals)
{
    std::cout << "Constraint: " << this->getName() << std::endl;
    std::cout << "Lower Bound: " << this->getLB() << ", Upper Bound: " << this->getUB() << std::endl;
    double eval = 0.0;
    for (size_t k = 0; k < Vars.size(); ++k) {
        NumVar* var = Vars[k];
		double val = Vals[k];
        double coef = this->getCoefficient(var);
		if (fabs(coef) > EPSILON && val > EPSILON) {
			std::cout << coef << var->getName() << "(" << val << ") +"; // Print coefficient, variable name, and value
        }
		eval += coef * val;
	}
	std::cout <<" = " << eval << std::endl;
}


 Env::~Env() {

     this->end(); // Clean up the CPLEX environment

     for (auto& c : _constraints) {
		 c->setasended();
         delete c;
     }
	 _constraints.clear();

     for (auto& v : _X) {
		 v->setasended();
         delete v;
     }
     _X.clear();

	 

 }

 void Env::addConstraint(Constraint* constraint)
 {
     constraint->toStandarForm();

     _model.add(*constraint);
     _constraints.push_back(constraint);

     // Every 1000 adds, restructure _constraints to avoid memory fragmentation
     if (_constraints.size() % 1000 == 0) {
         std::vector<Constraint*> tmp;
         tmp.reserve(_constraints.size());
         for (auto* c : _constraints) tmp.push_back(c);
         _constraints.swap(tmp);
     }
 }
  
 void Env::addConstraints(vector<Constraint*>& constraints)
 {
     for (auto c : constraints) addConstraint(c);
 }

 void Env::setupConstraint(Constraint* constraint) { constraint->addVars(_X); }

 void Env::removeConstraints(vector<Constraint*>& constraints)
 {
     // Step 1: Remove constraints from the model and _constraints vector
     for (Constraint* c : constraints) {
         // Remove from model
        try {
            _model.remove(*c);
        }
        catch (const IloException& e) {
            std::cerr << "Error removing constraint from model: " << e.getMessage() << std::endl;
        }
         // Remove from _constraints vector
         auto it = std::find(_constraints.begin(), _constraints.end(), c);
         if (it != _constraints.end()) {
             _constraints.erase(it);
         }

         // Step 2: Remove and delete slack variable if present
         NumVar* SlackVariable = c->getSlackVariable();
         if (SlackVariable != nullptr) {
             // Remove from model
             try {
                 _model.remove(*SlackVariable);
             }
             catch (const IloException& e) {
                 std::cerr << "Error removing slack variable from model: " << e.getMessage() << std::endl;
			 }

             // Remove from _X vector
             auto vit = std::find(_X.begin(), _X.end(), SlackVariable);
             if (vit != _X.end()) {
                 _X.erase(vit);
             }

             // Remove from _varIndexMap
             auto mapIt = _varIndexMap.find({SlackVariable->geti(), SlackVariable->getj(), SlackVariable->getType()});
             if (mapIt != _varIndexMap.end()) {
                 _varIndexMap.erase(mapIt);
             }

             // Delete
             delete SlackVariable;
             SlackVariable = nullptr;
         }

         // Delete constraint
         delete c;
     }
     // Clear the input vector to avoid dangling pointers
     constraints.clear();
 }

 void Env::addVar(NumVar* X) { 
     _X.push_back(X); 
     _varIndexMap[{X->geti(), X->getj(), X->getType()}] = _X.size() - 1; // Store the index of the variable in the map

     // Every 1000 adds, restructure _X to avoid memory fragmentation
     if (_X.size() % 1000 == 0) {
         // Move elements to a new vector to force memory compaction
         std::vector<NumVar*> tmp;
         tmp.reserve(_X.size());
         for (auto* v : _X) tmp.push_back(v);
         _X.swap(tmp);
         // No need to update _varIndexMap as indices remain the same
     }
 }

 NumVar* Env::getVar(int i, int j, Solver::VarType type) const {
     VarKey key = { i, j, type };
     auto it = _varIndexMap.find(key);
     if (it != _varIndexMap.end()) {
         return _X[it->second];
     }
     return nullptr; // Not found
 }


 double Env::getVarVal(NumVar& X) const { return _cplex.getValue(X); }



 IloNumArray Env::getVarVals(std::vector<NumVar*>& X)
 {
     IloNumVarArray vars(static_cast<IloEnv&>(*this));
     IloNumArray vals(static_cast<IloEnv&>(*this));

     for (auto x : X) vars.add(*x);
     try {
         _cplex.getValues(vals, vars);
     }
     catch (IloCplex::Exception e) {
         //cout << e.getMessage() << endl;
         vars.end();
         vals.end();
         throw e;
     }

     vars.end();
     return vals;
 }


 IloNumArray Env::getDuals()
 {
     // TODO: insert return statement here
     IloRangeArray ranges(*this);
     IloNumArray duals(*this);
     for (int k = 0; k < _constraints.size(); k++) ranges.add(*_constraints[k]);
     _cplex.getDuals(duals, ranges);
     ranges.end();

     return duals;
 }


 //inline bool Env::isSlackVar(int k) { return _X[k]->getType() == Solver::SlackVar; }

 //inline bool Env::isEdgeVar(int k) { return _X[k]->getType() == Solver::EdgeVar; }

 vector<double> Env::getInvRow(int k)
 {
     // Get the CPLEX environment and problem pointer
     CPXENVptr envPtr = _cplex.getImpl()->getCplexEnv();
     CPXLPptr lpPtr = _cplex.getImpl()->getCplexLp();

     std::vector<double> row(_cplex.getNcols());
     assert(_cplex.getNcols() == _X.size());

     // Call CPXbinvrow to get the k-th row of the basis inverse
     int status = CPXXbinvarow(envPtr, lpPtr, k, row.data());
     assert(status == 0);


     return row;
 }

 vector<double> Env::getInvBasisRow(int k)
 {
     // Get the CPLEX environment and problem pointer
     CPXENVptr envPtr = _cplex.getImpl()->getCplexEnv();
     CPXLPptr lpPtr = _cplex.getImpl()->getCplexLp();

     int numRows = _cplex.getNrows();
     std::vector<double> row(numRows);
     std::vector<double> rhs(numRows);
     if (numRows != _constraints.size())
         throw runtime_error("Invalid number of rows");

     // Call CPXbinvrow to get the k-th row of the basis inverse
     int status = CPXXbinvrow(envPtr, lpPtr, k, row.data());
     assert(status == 0);
     CPXXgetrhs(envPtr, lpPtr, rhs.data(), 0, numRows - 1);

     return row;
 }

 double Env::getInvRHS(int k)
 {
     // Get the CPLEX environment and problem pointer
     CPXENVptr envPtr = _cplex.getImpl()->getCplexEnv();
     CPXLPptr lpPtr = _cplex.getImpl()->getCplexLp();

     int numRows = _cplex.getNrows();
     std::vector<double> row(numRows);
     std::vector<double> rhs(numRows);
     if (numRows != _constraints.size())
         throw runtime_error("Invalid number of rows");

     // Call CPXbinvrow to get the k-th row of the basis inverse
     int status = CPXXbinvrow(envPtr, lpPtr, k, row.data());
     assert(status == 0);
     CPXXgetrhs(envPtr, lpPtr, rhs.data(), 0, numRows - 1);

     return ORUtils::dotproduct(row, rhs);
 }

 pair<int*, double*> Env::getBHead()
 {
     // Get the CPLEX environment and problem pointer
     CPXENVptr envPtr = _cplex.getImpl()->getCplexEnv();
     CPXLPptr lpPtr = _cplex.getImpl()->getCplexLp();

     int* head = new int[_cplex.getNrows()];
     double* x = new double[_cplex.getNrows()];

     CPXXgetbhead(envPtr, lpPtr, head, x);

     return make_pair(head, x);
 }

 IloCplex::BasisStatusArray Env::getBasisStatus(const std::vector<NumVar*>& X)
 {
     IloCplex::BasisStatusArray cstat(*this);
     IloNumVarArray Var(*this);
     for (auto x : X) Var.add(*x);

     _cplex.getBasisStatuses(cstat, Var);
     Var.end();
     return cstat;
 }

 std::vector<std::string> Env::getSortedColNames() {
     // Get the CPLEX environment and problem pointer
     CPXENVptr envPtr = _cplex.getImpl()->getCplexEnv();
     CPXLPptr lpPtr = _cplex.getImpl()->getCplexLp();

     int ncols = CPXXgetnumcols(envPtr, lpPtr);
     if (ncols <= 0) return {};

     // Prepare storage for names
     CPXSIZE surplus = 0;
     CPXSIZE namebufsize = 0;
     int status = CPXXgetcolname(envPtr, lpPtr, NULL, NULL, 0, &namebufsize, 0, ncols - 1);
     //if (status != 0 || namebufsize >= 0) return {}; // error or no names

     std::vector<char> namebuf(static_cast<size_t>(-namebufsize));
     std::vector<char*> nameptrs(ncols, nullptr);
     status = CPXXgetcolname(envPtr, lpPtr, nameptrs.data(), namebuf.data(), -namebufsize, &surplus, 0, ncols - 1);
     if (status != 0) return {};

     std::vector<std::string> names;
     for (int i = 0; i < ncols; ++i) {
         if (nameptrs[i])
             names.emplace_back(nameptrs[i]);
         else
             names.emplace_back(""); // fallback for unnamed
     }
     return names;
 }

 //We compute the reduced cost of a variable

 inline double Env::getReducedCost(NumVar* x) { 
	 try {
		 return _cplex.getReducedCost(*x);
	 }
	 catch (IloCplex::Exception e) {
		 throw e;
	 }
 }

 void Env::addCol(NumColumn* col)
 {
     int nextid = int(_X.size());
     IloNumColumn iloCol(*this);
     iloCol += _obj(col->getObj());
     //for (int k = 0; k < _constraints.size(); k++) if(fabs(C[k])>CGEnv::EPSILON) iloCol += (*_constraints[k])(C[k]);
     for (auto& it : *col) iloCol += (*_constraints[it.first])(it.second);
     NumVar* x = new NumVar(*this, iloCol, col->geti(), col->getj());
     assert(nextid == x->getId());
     _X.push_back(x);
     _model.add(*x);

     col->setInModel(x);
     col->clear();

     iloCol.end();
 }

 double Env::getDistanceToCut(Constraint& R)
 {

     //IloNumArray Vals = getVarVals(); //We need to compute the values of the variables first
     //double val = 0.0;
     //for (int i = 0; i < _X.size(); i++) {
     //    if (Vals[i] > EPSILON) {
     //        double coef = R.getCoefficient(_X[i]);
     //        val += Vals[i] * coef;
     //    }
     //}
	 
     
     IloExpr expr(R.getExpr());
     IloNum lhs = R.getLB();
     IloNum rhs = R.getUB();

     IloNum val = _cplex.getValue(expr);

     double diff = fmax((lhs - val), (val - rhs));
     //if (diff < 0) cout <<"Csts: " << expr << endl;
     //assert(diff >= 0.0);
     expr.end();
     //Vals.end();
     return diff;

     //IloExpr::LinearIterator it = expr.getLinearIterator();
     //double qsum = 0.0;
     //while (it.ok()) {
     //    if (it.getVar().getUB() > CGEnv::EPSILON)
     //        qsum += pow(it.getCoef(), 2);
     //    ++it;
     //}
     //if (qsum < CGEnv::EPSILON) {
     //    cout << R << endl;
     //}
     //if (qsum < CGEnv::EPSILON) return NULL;
     //assert(qsum > CGEnv::EPSILON);
     //qsum = sqrt(qsum);

     //double dist = diff / qsum;
     ////cout << _R->getName() <<  ": " << diff << " - " << qsum << " - " << dist << endl;
     ////assert(fabs(diff) >= fabs(dist));

     //expr.end();
     //return dist;
 }

 void Env::setObjVarCoef(NumVar* X, double obj) { _obj.setLinearCoef(*X, obj); }

 void Env::setObjective(std::function<double(NumVar*)> costfunction)
 {
     IloExpr expr(*this);
     for (auto& x : _X) {
         double coef = costfunction(x);
         if (fabs(coef) > EPSILON) expr += (*x) * coef;
     }
     _obj.setExpr(expr);
     expr.end();
 }

 double Env::getGap() {
     double g = 1.0;
     try {
         g = _cplex.getMIPRelativeGap();
     }
     catch (IloException& e) {
         cout << "Error: " << e.getMessage() << endl;
     }
     return g;
 }

 Solver::SolverStatus Env::runsolver()
 {
     
     if (_cplex.solve()) {
		 auto status = _cplex.getStatus();
		 //cout << "Solution status = " << _cplex.getStatus() << endl;
		 //cout << "Solution value  = " << _cplex.getObjValue() << endl;
		 if (status == IloAlgorithm::Optimal) return Solver::SolverStatus::Optimal;
		 else if (status == IloAlgorithm::Infeasible) return Solver::SolverStatus::Infeasible;
		 else return Solver::SolverStatus::Unknown;
     }
	 else return Solver::SolverStatus::Infeasible;
 }

 Solution* Env::recoverSolution()
 {

     IloNumVarArray vars(static_cast<IloEnv&>(*this));
     IloNumArray vals(static_cast<IloEnv&>(*this));

     for (auto x : _X) if (x->getType() == Solver::EdgeVar) vars.add(*x);
     try {
         _cplex.getValues(vals, vars);
     }
     catch (IloCplex::Exception e) {
         cout << e.getMessage() << endl;
         vars.end();
         vals.end();
		 return NULL;
     }

     Solution*  S = new Solution(_I);
    
     for (int k = 0; k < vars.getSize(); k++) {
         if (vals[k] > 0) {
             int i = _X[k]->geti();
             int j = _X[k]->getj();
             double val = vals[k];
             S->addEdge(i, j, val);
         }
     }
     vars.end();
     vals.end();
     return S;

     

  //   Solution* S = NULL;

	 //try {
  //       S = new Solution(_I);
		// IloNumArray vals = getEdgeVarVals();

		// for (int k = 0; k < _X.size(); k++) {
		//	 if (/*isEdgeVar(k) &&*/ vals[k]>EPSILON) {
		//		 int i = _X[k]->geti();
		//		 int j = _X[k]->getj();
		//		 double val = vals[k];
		//		 S->addEdge(i, j,val);
		//	 }
		// }
		// vals.end();
  //       return S;
	 //}
	 //catch (exception& e) {
		// cout << "Error: " << e.what() << endl;
		// if (S != NULL) delete S;
  //       S = NULL;
	 //}
  //   catch (IloCplex::Exception e) {
  //       cout << e.getMessage() << endl;
  //       if (S != NULL) delete S;
  //       S = NULL;
  //   }

  //   return S;
 }

 void Env::convertToIP()
 {
	 IloNumVarArray vars(*this);
	 for (auto& x : _X) {
		 if (x->getType() == Solver::EdgeVar)
			 vars.add(*x);
	 }

	 IloConversion c(*this, vars, ILOINT);
	 _model.add(c);
 }

 void Env::convertToLP()
 {
	 IloNumVarArray vars(*this);
	 for (auto& x : _X) {
		 if (x->getType() == Solver::EdgeVar)
			 vars.add(*x);
	 }
	 IloConversion c(*this, vars, ILOFLOAT);
	 _model.add(c);
 }

 void Env::DebugCheckConstraints(const std::vector<Constraint*>& C)
 {
	 string solfile = _I->knownsolpath();

     Solution* S = NULL;

     try
     {
		 S = new Solution(solfile, _I);
     }
     catch (const std::exception& e)
     {
		 cout << "Error reading solution file: " << e.what() << endl;
		 if (S != NULL) delete S;
		 S = NULL;
		 return;
     }

	 auto Vals = extractVarVals(S);
	 auto& Vars = _X; // Get the variables from the environment

	 //We check if all constraints are satisfied
	 bool allSatisfied = true;
     for (auto& c : C) {
         if (!c->isValidUnderConstraint(_X,Vals)) {
             cout << "Constraint not satisfied: " << c->getName() << endl;
			 c->print(Vars, Vals);
			 allSatisfied = false;
         }
	 }

     if (!allSatisfied) {
         exportModel("DebugModel.lp");
		 throw runtime_error("Some constraints are not satisfied by the solution.");
     }
     else
		 cout << "All constraints are satisfied by the best known solution." << endl;


     delete S;

 }

 std::vector<double> Env::extractVarVals(Solution* S) const
 {
     // Pseudocode:
     // 1. Map Solution's flow edges to EdgeVar variables in _X.
     // 2. For each edge in S, set the corresponding EdgeVar value in the output array.
     // 3. For CapVar variables, propagate from depot (assume node 0) with initial capvar=0, 
     //    and for each outgoing edge, set capvar at target as capvar at source + flow.
     // 4. Fill the output array for CapVar variables accordingly.
     // 5. All other variables set to 0.

     std::vector<double> vals(_X.size(), 0.0);

     // 1. Map edge flows
     std::map<std::pair<int, int>, double> edgeFlow;
     for (auto& e : *S) {
         edgeFlow[{e->geti(), e->getj()}] = e->getFlow();
     }

     // 2. Set EdgeVar values
     for (size_t k = 0; k < _X.size(); ++k) {
         NumVar* var = _X[k];
         if (var->getType() == Solver::EdgeVar) {
             auto it = edgeFlow.find({var->geti(), var->getj()});
             if (it != edgeFlow.end()) {
                 vals[k] = it->second;
             }
         }
     }

     // 3. CapVar propagation
     // Find all CapVar variables and build adjacency
     std::map<int, std::vector<std::pair<int, double>>> outEdges; // node -> list of (to, flow)
     for (auto& ef : edgeFlow) {
         int from = ef.first.first;
         int to = ef.first.second;
         double flow = ef.second;
         outEdges[from].emplace_back(to, flow);
     }

     // Find all CapVar indices
     std::map<int, size_t> capVarIdx; // node -> index in _X
     for (size_t k = 0; k < _X.size(); ++k) {
         NumVar* var = _X[k];
         if (var->getType() == Solver::CapVar) {
             capVarIdx[var->geti()] = k;
         }
     }

     // Forward propagate from depot (assume node 0)
     std::map<int, double> capVal;
     std::queue<int> q;
     capVal[0] = 0.0;
     q.push(0);
     while (!q.empty()) {
         int u = q.front(); q.pop();
         double cu = capVal[u];
         // Set CapVar value if exists
         auto it = capVarIdx.find(u);
         if (it != capVarIdx.end()) {
             vals[it->second] = cu;
         }
         for (auto& [v, f] : outEdges[u]) {
             if (capVal.find(v) == capVal.end() && f>0) {
                 // Add demand of node u when propagating
                 double demand_v = _I->Node(v)->demand();
                 capVal[v] = cu + demand_v;
                 q.push(v);
             }
         }
     }

     // Clear unused mapping
	 outEdges.clear();
	 capVal.clear();
	 capVarIdx.clear();

     // 5. Eval slack variables in each constraints
     for(auto c : _constraints) {
         NumVar* slackVar = c->getSlackVariable();
         if (slackVar != NULL) {
             double coef = c->getSlackCoefficient();
             double eval = 0.0;
             for(size_t k = 0; k < _X.size(); ++k) {
                 NumVar* var = _X[k];
                 double varcoef = c->getCoefficient(var);
                 eval += varcoef * vals[k]; // Use vals[k] for the current value
			 }
             double slackVal = 0.0;
             if(coef > 0.5) { // Positive slack variable
				 slackVal = c->getUB() - eval; // slack = UB - eval
             }
             else { // Negative slack variable
				 slackVal = eval - c->getLB(); // slack = eval - LB
			 }

             if (slackVal < 0) {
                 c->print(_X, vals);
				 throw runtime_error("Negative slack value detected: " + to_string(slackVal) + " for constraint: " + c->getName());
             }

			 //We get the index of the slack variable in _X
			 auto it = _varIndexMap.find({ slackVar->geti(), slackVar->getj(), slackVar->getType() });
             if (it != _varIndexMap.end()) {
                 int _varIndex = it->second;
				 // Set the value in vals
				 vals[_varIndex] = slackVal;
             }
             else 
				 throw runtime_error("Slack variable not found in _varIndexMap");

         }
	 }




     return vals;
 }
