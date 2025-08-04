#include "Solution.hpp"



Solution::Solution(Instance* I) : vector<FlowEdge*>(), _I(I) {
    _pred.resize(_I->nnodes(), -1);
    _adj.resize(_I->nnodes());
}

Solution::Solution(string filename, Instance* I) : vector<FlowEdge*>(), _I(I) {

    _pred.resize(_I->nnodes(), -1);
    _adj.resize(_I->nnodes());
   
    ifstream file(filename);

    if (!file.is_open()) {
        throw invalid_argument("Error: could not open file " + filename);
    }

    int k = 1;
    while (ORUtils::gotoKeyword(file, "#" + to_string(k) + ":")) {
        int i = _I->InitialDepotId();
        int j = i;
        string line;
        getline(file, line);
        istringstream iss(line);
        while (iss >> j) {
			addEdge(i, j);
            i = j;
        }
		addEdge(j, _I->FinalDepotId());
        k++;
    }

    if (ORUtils::gotoKeyword(file, "Cost")) {
        double cost;
        file >> cost;
        //cout << "Cost: " << cost << endl;
    }

    file.close();
}

Solution::Solution(Solution* S)
{
	_I = S->_I;
	_pred = S->_pred;
	_adj = S->_adj;
	for (auto e : *S) {
		FlowEdge* e2 = new FlowEdge();
		e2->id = e->id;
		e2->i = e->i;
		e2->j = e->j;
		e2->rf = e->rf;
		e2->f = e->f;
		e2->r = e->r;
		this->push_back(e2);
	}
}

Solution::~Solution()
{
	for (auto e : *this) delete e;
}

int Solution::addEdge(int i, int j, double fij, double fji) {
	//We check first id the edge already exists
	for (auto k : _adj[i]) {
		if (this->at(k)->j == j) {
			this->at(k)->f += fij;
			return this->at(k)->id;
		}
	}


    FlowEdge* eij = new FlowEdge();
    FlowEdge* eji = new FlowEdge();
    eij->id = (int)this->size();
    eij->i = i;
    eij->j = j;
    eij->f = fij;
    eij->rf = 0;
    eji->id = (int)this->size() + 1;
    eji->i = j;
    eji->j = i;
    eji->f = fji;
    eji->rf = 0;
    eij->r = eji;
    eji->r = eij;
    this->push_back(eij);
    this->push_back(eji);
    _adj[i].push_back(eij->id);
    _adj[j].push_back(eji->id);

	return eij->id;
}

double Solution::getFlow(int i, int j)
{
	if (_adj[i].empty()) return 0.0;
	auto it = find_if(_adj[i].begin(), _adj[i].end(), [j, this](int k) {return this->at(k)->j == j; });
	if (it != _adj[i].end()) return this->at(*it)->f;
	else return 0.0;
}

double Solution::getFlow(const std::vector<int>& sources, const std::vector<int>& sinks)
{
    std::vector<bool> inSource(_I->nnodes(), false);
    for (int i : sources) inSource[i] = true;
    std::vector<bool> inSink(_I->nnodes(), false);
    for (int i : sinks) inSink[i] = true;

    double flow = 0.0;
    for (FlowEdge* e : *this) {
        if (inSource[e->geti()] && inSink[e->getj()])
            flow += e->f;
    }
    

        
    return flow;
}

vector<FlowEdge*> Solution::getAdjacentEdges(int i) const
{
	if (_adj[i].empty()) return vector<FlowEdge*>();
	vector<FlowEdge*> edges;
	for (auto k : _adj[i]) {
		edges.push_back(this->at(k));
	}
	return edges;
}

double Solution::fitness() {
    if (this->empty()) return 0;

	//auto R = extractroutes();
	//double cost = 0;
	//for (auto r : R) {
	//	for (int i = 0; i < r.size() - 1; i++) {
	//		cost += _I->cost(r[i], r[i + 1]);
	//	}
	//}
	//return cost;

    return accumulate(begin(), end(), 0.0, [&](double sum, FlowEdge* e) {return sum + e->f * _I->cost(e->i, e->j); });
}

bool Solution::isInteger() {
    for (auto e : *this) {
        if (!ORUtils::isInteger(e->f)) return false;
        //if ((->rf > EPSILON || e->rf <= 1.0 - EPSILON) return false;
    }
    return true;
}

vector<FlowEdge*> Solution::getPath(const vector<int>& sources, const vector<int>& sinks)
{
    vector<FlowEdge*> Path;

    fill(_pred.begin(), _pred.end(), -2);

    for (auto t : sinks) _pred[t] = -3;

    vector<int> Q;
    for (auto s : sources) {
        Q.push_back(s);
        _pred[s] = -1;
    }


    while (!Q.empty()) {
        int i = Q.back();
        Q.pop_back();

        for (int k : _adj[i]) {
			FlowEdge* e = this->at(k);
            int j = e->j;

            if (e->rf < EPSILON) continue;

            //We check if we reach one of the targets
            if (_pred[j] == -3) {
                _pred[j] = e->id;
                Q.clear();
                for (int k = _pred[j]; k >= 0; k = _pred[this->at(k)->i]) Path.push_back(this->at(k));
                reverse(Path.begin(), Path.end());
                break;
            }
            else if (_pred[j] == -2) {
                _pred[j] = e->id;
                Q.push_back(j);
            }
        }
    }


    return Path;
}

void Solution::resetResidualFlow() {
    for (auto e : *this) e->rf = e->f;
    fill(_pred.begin(), _pred.end(), -1);
}

std::vector<int> Solution::reachableNodes(const std::vector<int>& sources)
{

	std::vector<bool> visited(_I->nnodes(), false);
	for (int s : sources) visited[s] = true;

	std::vector<int> queue(sources.begin(), sources.end());
	while (!queue.empty()) {
		int u = queue.back();
		queue.pop_back();
		for (int k : _adj[u]) {
			FlowEdge* e = this->at(k);
			int v = e->j;
			if (e->rf > EPSILON && !visited[v]) {
				visited[v] = true;
				queue.push_back(v);
			}
		}
	}

	std::vector<int> components;
	for (int i = 0; i < _I->nnodes(); ++i) {
		if (visited[i]) components.push_back(i);
	}

	return components; // All nodes are connected
}

vector<vector<int>> Solution::extractroutes()
{
	
    vector<int> firstnodes;

	int origin = _I->InitialDepotId();
	int destination = _I->FinalDepotId();
	
	resetResidualFlow();

	vector<vector<int>> routes;
    vector<FlowEdge*> Path = getPath({ origin }, { destination });
    while (!Path.empty()) {
		
		for (auto e : Path) {
			e->rf -= 1;
			e->r->rf += 1;
		}

		vector<int> route;
		for (auto e : Path) {
			route.push_back(e->i);
		}
		route.push_back(Path.back()->j);
		routes.push_back(route);
		
        Path = getPath({ origin }, { destination });
	}

    return routes;
}

bool Solution::isFeasible(bool displayerror) {
    vector<string> Error;

    if (!isInteger()) {
        Error.push_back("Solution is not integer");
    }
    else {
        auto routes = extractroutes();
        if (routes.size() > _I->FleetSize()) {
            Error.push_back("Number of routes is greater than the fleet size");
        }
        for (auto& r : routes) {
            double Q = 0.0;
            for (auto i : r) {
                Q += _I->Node(i)->demand();
            }

            if (Q > _I->VehicleCap()) {
                Error.push_back("Route with demand greater than the vehicle capacity");
            }
        }
    }

    vector<double> visited(_I->nnodes(), 0);
	for (auto e : *this) {	
	    visited[e->i] += e->f/2.0;
        visited[e->j] += e->f/2.0;
	}

    if (visited[_I->FinalDepotId()] > _I->FleetSize()) {
		Error.push_back("Final depot is visited more than the fleet size");
	}
    visited[_I->FinalDepotId()] = 1;
	visited[_I->InitialDepotId()] = 1;

    double nvisits = accumulate(visited.begin(), visited.end(), 0.0, [](double sum, double v) -> double {return sum + v; });
    
    //cout << fabs(nvisits - (double)_I->nnodes()) << endl;

    if (fabs(nvisits-(double)_I->nnodes())>EPSILON ) {
		//exportGraph("solution_graph.dot");
        Error.push_back("Number of nodes is different from the number of nodes in the instance");
    }
    for (int i = 0; i < _I->nnodes(); i++) {
        if (visited[i] < EPSILON) {
            Error.push_back("Node " + to_string(i) + " is missing");
        }
		else if (visited[i] > 1.0 + EPSILON) {
            Error.push_back("Node " + to_string(i) + " is visited more than once (" + to_string(visited[i]) + ")");
		}
        else if (visited[i] < 1.0 - EPSILON) {
			Error.push_back("Node " + to_string(i) + " is visited partially");
		}
    }
    if (Error.size() > 0) {
        if (displayerror) {
            cout << "Solution is not feasible" << endl;
            for (auto& e : Error) {
                cout << e << endl;
            }
        }
        return false;
    }
    else return true;
}

void Solution::print() {
    cout << "Solution: ";
	bool isint = isInteger();
	if (isint) cout << "Integer solution" << endl;
	else cout << "Fractional solution" << endl;

    //isint = false;
    if (isint) {
        auto routes = extractroutes();
        reverse(routes.begin(), routes.end());
        int k = 1;
        for (auto& r : routes) {
            cout << "\t Route #" << k++ << ": ";
            for (auto i : r) {
                cout << _I->Node(i)->id() << " ";
            }
            cout << endl;
        }
    }
    else {
        for (auto e : *this) {
            if (e->f > EPSILON) {
                cout << "\t" << e->i << " -> " << e->j << " : " << e->f << endl;
            }
        }
    }
    
    cout << "Cost " << fitness() << endl;
    cout << endl;
}

void Solution::exportGraph(const std::string& filename)
{
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "Error: could not open file " << filename << " for writing." << std::endl;
        return;
    }

    ofs << "digraph SolutionGraph {\n";
    ofs << "    node [shape=circle];\n";

    // Output nodes with labels and distinctive color for depots
    for (int i = 0; i < _I->nnodes(); ++i) {
        const node* n = _I->Node(i);
        ofs << "    " << i << " [label=\"" << n->id();
        if (n->isInitialDepot()) ofs << " (Start)";
        else if (n->isFinalDepot()) ofs << " (End)";
        ofs << "\"";
        if (n->isInitialDepot()) {
            ofs << ", style=filled, fillcolor=lightgreen";
        }
        else if (n->isFinalDepot()) {
            ofs << ", style=filled, fillcolor=lightblue";
        }
        ofs << "];\n";
    }

    // Output edges with flow as label
    for (const FlowEdge* e : *this) {
        if (e->f > EPSILON) {
            ofs << "    " << e->i << " -> " << e->j << " [label=\"" << e->f << "\"";
            if (e->f > 0.0) ofs << ", penwidth=2";
            ofs << "];\n";
        }
    }

    ofs << "}\n";
    ofs.close();
}

Solution* Solution::readSolution(const std::string& filename, Instance* I)
{
	Solution* sol = nullptr;
    try
    {
		sol = new Solution(filename, I);
    }
    catch (const std::exception& e)
    {
		std::cout << e.what() << std::endl;
		if (sol != nullptr) delete sol;
		sol = nullptr;
    }
    return sol;
}

std::pair<std::vector<int>, std::vector<int>> Solution::maxFlow(const std::vector<int>& sources, const std::vector<int>& sinks)
{
    // Pseudocode:
    // 1. Reset all residual flows to the original flow values.
    // 2. While there exists an augmenting path from sources to sinks (using getPath):
    //    a. Find the minimum residual flow along the path (bottleneck).
    //    b. Augment the flow along the path by the bottleneck value.
    //    c. Update the residual flows accordingly.
    // 3. After no more augmenting paths, collect the set of reachable nodes from sources in the residual graph (S).
    // 4. The rest are in T.
    // 5. Return the pair (S, T).

    resetResidualFlow();

    // Step 2: Edmonds-Karp style max-flow
	double flow = 0.0;
    while (true) {
        std::vector<FlowEdge*> path = getPath(sources, sinks);
        if (path.empty()) break;

        // Find bottleneck
        double bottleneck = std::numeric_limits<double>::max();
        for (auto e : path) {
            if (e->rf < bottleneck) bottleneck = e->rf;
        }

		assert(bottleneck > EPSILON && "Bottleneck should be greater than EPSILON.");

        // Augment flow along the path
        for (auto e : path) {
            e->rf -= bottleneck;
            e->r->rf += bottleneck;
        }

		flow += bottleneck;
    }


    // Step 3: BFS from sources in residual graph to find reachable nodes (S)
    std::vector<bool> visited(_I->nnodes(), false);
    std::vector<int> queue;
    for (int s : sources) {
        queue.push_back(s);
        visited[s] = true;
    }
    while (!queue.empty()){
		int u = queue.back();
		queue.pop_back();
        for (int k : _adj[u]) {
            FlowEdge* e = this->at(k);
            int v = e->j;
            if (e->rf > EPSILON && !visited[v]) {
                visited[v] = true;
                queue.push_back(v);
            }
        }
    }

    std::vector<int> S, T;
    for (int i = 0; i < _I->nnodes(); ++i) {
        if (visited[i]) S.push_back(i);
        else T.push_back(i);
    }

    if(fabs(flow - getFlow(S, T)) >= EPSILON)
		exportGraph("max_flow_graph.dot");

	assert(fabs(flow - getFlow(S, T)) < EPSILON && "Max flow does not match the computed flow.");

    return std::make_pair(S, T);
}

