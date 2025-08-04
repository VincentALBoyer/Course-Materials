#include "IPSolver.hpp"
#include "../Solution.hpp"
#include <queue>
#include <set>
#include <functional>

unsigned int BBNode::_index = 0;



void BBNode::applyVarStatus(NumVar* X, Solver::VarStatus status) {
#if DEBUG
	if (X->getType() != Env::VarType::EdgeVar) throw std::exception("Variable is not an edge variable");
#endif
	X->setLB(0.0);
	X->setUB(1.0);
	if (status == Solver::VarStatus::Down) X->setUB(0.0);
	else if (status == Solver::VarStatus::Up) X->setLB(1.0);
}

void BBNode::revokeVarStatus(NumVar* X)
{
#if DEBUG
	if (X->getType() != Env::VarType::EdgeVar) throw std::exception("Variable is not an edge variable");
#endif
	X->setLB(0.0);
	X->setUB(1.0);
}

void BBNode::applyNodeVarStatus(BBNode* node)
{
	BBNode* parent = node;
	while (parent != NULL) {
		for (auto& it : parent->_statusMap) {
			NumVar* X = it.first;
			Solver::VarStatus status = it.second;
			applyVarStatus(X, status);
		}
		parent = parent->_parent;
	}
}

void BBNode::revokeNodeVarStatus(BBNode* node)
{
	BBNode* parent = node;
	while (parent != NULL) {
		for (auto& it : parent->_statusMap) {
			NumVar* X = it.first;
			revokeVarStatus(X);
		}
		parent = parent->_parent;
	}
}

//Constructor for the Root Node


BBNode::~BBNode() {
	//Delete Children
	if (_left != NULL)
		delete _left;
	if (_right != NULL)
		delete _right;

	//Delete Parent
	if (_parent != NULL) {
		if (_parent->_left == this) _parent->_left = NULL;
		if (_parent->_right == this) _parent->_right = NULL;
		if (_parent->_left == NULL && _parent->_right == NULL) delete _parent;
	}

}

void BBNode::generateChildren(NumVar* X) {
	_left = new BBNode(X, Solver::VarStatus::Down, this, _objval);
	_right = new BBNode(X, Solver::VarStatus::Up, this, _objval);
}

void BBNode::generateChildren(VarStatusMap& statusMapLeft, VarStatusMap& statusMapRight)
{
    _left = new BBNode(statusMapLeft, this, _objval);
	_right = new BBNode(statusMapRight, this, _objval);
}

void BBNode::setStatus(NumVar* X, Solver::VarStatus status) {
#if DEBUG
	if (_statusMap.find(X) != _statusMap.end()) throw std::exception("Variable already set");
#endif
	_statusMap[X] = status;
}

void BBNode::setLocalVarBound() {
	for (auto& it : _statusMap) {
		Solver::VarStatus status = it.second;
		NumVar* X = it.first;
		if (status == Solver::VarStatus::Down) X->setLB(0.0);
		else if (status == Solver::VarStatus::Up) X->setUB(0.0);
	}
}

void BBSolver::solvemethod(Solution* S)
{
    const double rounding_fraction = 0.9; // 15% of total variables as threshold
    if (S != NULL) {
        updateBestTour(S);
    }

    _LP->setparam(Solver::TimLim, timeLeft());

    if (_BBNodes.empty()) {
        //We create the root node
        BBNode* root = new BBNode();
        addNode(root);
    }

    if (!_BBNodes.empty())
        printlog(0, _BBNodes.top());

    int iter = 0;
    while (!_BBNodes.empty() && gap() > EPSILON && timeLeft() > EPSILON) {
        iter++;

        //We get the node with the best objective value
        BBNode* node = _BBNodes.top();
        _BBNodes.pop();

        //We apply the variable status
        BBNode::applyNodeVarStatus(node);
        //We solve the LP
        _LP->solve();

        //We update the Node Status
        auto status = _LP->getStatus();
        if (status == Solver::SolverStatus::Infeasible) {
            node->setStatus(Solver::NodeStatus::InfeasibleNode);
        }
        else if (status == Solver::SolverStatus::Optimal) {
            node->setStatus(Solver::NodeStatus::FeasibleNode);

            //We get the solution
            Solution* NodeSol = _LP->recoversolution();

            //We update the objective value
            node->setObjVal(NodeSol->fitness());

            //We check if the solution is integer
            vector<FlowEdge*> candidates = getBranchingCandidates(NodeSol);
            double frac_percent = (NodeSol->size() > 0) ? (double)candidates.size() / ((double)NodeSol->size() - (double)candidates.size()) : 0.0;
            if (node->getObjVal() >= _bestUB) {
                node->setStatus(Solver::NodeStatus::PrunedNode);
            }
            else if (NodeSol->isInteger()) {
                //We update the best solution
                updateBestTour(NodeSol);
                node->setStatus(Solver::NodeStatus::PrunedNode);
            }
            else{
                if(candidates.empty()) 
                    throw std::exception("No candidates for branching found");
                //cout << "Perc Frac: " << frac_percent << endl;
                if (frac_percent <= rounding_fraction) {
                    // Simple greedy heuristic
                    Solution* greedySol = greedyConstructSolution(NodeSol);
                    if (greedySol->isInteger() && greedySol->isFeasible(false) && greedySol->fitness() < _bestUB) {
                        updateBestTour(greedySol);
                    }
                    delete greedySol;
                }
                
                //FlowEdge* e = *min_element(candidates.begin(), candidates.end(), [](FlowEdge* e1, FlowEdge* e2) {return ORUtils::Fractionality(e1->getFlow()) < ORUtils::Fractionality(e2->getFlow()); });
				//FlowEdge* e = *max_element(candidates.begin(), candidates.end(), [this](FlowEdge* e1, FlowEdge* e2) {return _I->cost(e1->geti(), e1->getj()) < _I->cost(e2->geti(), e2->getj()); });
                FlowEdge* e = *min_element(candidates.begin(), candidates.end(), [this](FlowEdge* e1, FlowEdge* e2) {
					double f1 = ORUtils::Fractionality(e1->getFlow());
					double f2 = ORUtils::Fractionality(e2->getFlow());
					if (fabs(f1 - f2) < EPSILON) {
                        return _I->cost(e1->geti(), e1->getj()) < _I->cost(e2->geti(), e2->getj());
                    }
                    else return f1 < f2;
					
                    });

                
                //node->generateChildren(_LP->getVar(e));
                

				auto VarSplit = _splitVars(NodeSol, e->geti()); // Split the variables for the branching edge
				VarStatusMap statusMapLeft, statusMapRight;
                for (auto& var : VarSplit.first) {
                    if(var->getUB() >= EPSILON) // Only set status if the variable is not fixed
						statusMapLeft[var] = Solver::VarStatus::Down; // Left child takes Down status
                }
                for (auto& var : VarSplit.second) {
                    if (var->getUB() >= EPSILON) // Only set status if the variable is not fixed
					    statusMapRight[var] = Solver::VarStatus::Down; // Right child takes Down status on the complementary variables
				}
				node->generateChildren(statusMapLeft, statusMapRight);
				

                addNode(node->getLeft());
                addNode(node->getRight());
                
            }

            delete NodeSol;

        }
        else {
            throw std::exception("Unknown status");
        }
        
        //We revoke the variable status
        BBNode::revokeNodeVarStatus(node);

        //We print some logs: gap, bestobjval, bestbound
        printlog(1000, node);

        //We delete the node
        if(node->getStatus()== Solver::NodeStatus::InfeasibleNode || node->getStatus()== Solver::NodeStatus::PrunedNode ) delete node;
    }

    return;
}

void BBSolver::addNode(BBNode* node)
{
	// We set the pending status
	node->setStatus(Solver::NodeStatus::PendingNode);
	_BBNodes.push(node);
}

bool BBSolver::updateBestTour(Solution* S)
{
	if (!S->isInteger()) return false;
	if (_bestsol == NULL || S->fitness() < _bestsol->fitness()) {
		if (_bestsol != NULL) delete _bestsol;
		_bestsol = new Solution(S);
		_bestUB = _bestsol->fitness();
		return true;
	}
	else return false;
}

vector<FlowEdge*> BBSolver::getBranchingCandidates(Solution* S)
{
	//We get the fractional edges
	vector<FlowEdge*> candidates;
	for (auto e : *S) {
		if (ORUtils::isFractional(e->getFlow())) candidates.push_back(e);
	}
	return candidates;
}

void BBSolver::printlog(int freq, BBNode* N)
{
	static unsigned int iter = 0;

	int w = 20;
	if (freq <= 0) {
		cout << setw(w) << "Left" << setw(w) << "NodeId" << setw(w) << "Obj" << setw(w) << "UB" << setw(w) << "LB" << setw(w) << "Gap" << endl;
	}
	else {
		iter++;
		if (iter % freq == 0) {
			string UB = "n/a";
			string G = "n/a";
			if (_bestsol != NULL) {
				UB = to_string(_bestsol->fitness());
				G = to_string(gap());
			}
			if (!_BBNodes.empty()) {
				cout << setw(w) << _BBNodes.size() << setw(w) << N->getId() << setw(w) << N->getObjVal() << setw(w) << UB << setw(w) << _BBNodes.top()->getObjVal() << setw(w) << G << endl;
			}
		}
	}
}

Solution* BBSolver::recoversolution()
{
	if (_bestsol == NULL) return NULL;
	return new Solution(_bestsol);
}

// Greedy repair using NodeSol's flows
Solution* BBSolver::greedyConstructSolution(Solution* NodeSol) {
    int n = _I->nnodes();
    int K = _I->FleetSize();
    double Q = _I->VehicleCap();
    std::vector<bool> visited(n, false);
    visited[_I->InitialDepotId()] = true;
    visited[_I->FinalDepotId()] = true;
    
    // Store routes for each vehicle
    std::vector<std::vector<int>> routes(K);
    std::vector<double> route_caps(K, Q);
    double total_cost = 0.0;


    for (int k = 0; k < K; ++k) {
        int curr = _I->InitialDepotId();
        double cap = Q;
        std::vector<int>& route = routes[k];
        route.push_back(curr);
        while (true) {
            int next = -1;
            
            // Find the outgoing edge from curr with highest flow in NodeSol
            std::vector<FlowEdge*> edges = NodeSol->getAdjacentEdges(curr);
            // double best_flow = -1.0;
            //for (auto e : edges) {
            //    int j = e->getj();
            //    if (!visited[j] && j != _I->FinalDepotId()) {
            //        double demand = _I->Node(j)->demand();
            //        if (demand <= cap && e->getFlow() > best_flow) {
            //            best_flow = e->getFlow();
            //            next = j;

            //            if (best_flow >= 0.5) break;
            //        }
            //    }
            //}

            // Remove infeasible edges: visited, final depot, or demand > cap
            edges.erase(
                std::remove_if(edges.begin(), edges.end(), [&](FlowEdge* e) {
                    int j = e->getj();
                    double demand = _I->Node(j)->demand();
                    return visited[j] || j == _I->FinalDepotId() || demand > cap;
                }),
                edges.end()
            );
            if(!edges.empty()) {
                std::vector <double> flows;
                for (auto e : edges) flows.push_back(e->getFlow());
                size_t index = ORRandom::roulette<std::vector <double>>(flows);
				next = edges[index]->getj();
			}
			

            if (next == -1) break; // No more nodes can be added
            cap -= _I->Node(next)->demand();
            total_cost += _I->cost(curr, next);
            visited[next] = true;
            curr = next;
            route.push_back(curr);


        }
        // Go to final depot
        total_cost += _I->cost(curr, _I->FinalDepotId());
        route.push_back(_I->FinalDepotId());
        route_caps[k] = cap;
    }

    if(total_cost >= _bestUB) {
        // If the cost is already worse than the best known solution, return empty solution
        return new Solution(_I);
    }

    // Priority queue for unvisited nodes (higher demand first)
    using NodeDemand = std::pair<double, int>; // (demand, node)
    auto cmp = [](const NodeDemand& a, const NodeDemand& b) { return a.first < b.first; };
    std::priority_queue<NodeDemand, std::vector<NodeDemand>, decltype(cmp)> pq(cmp);
    for (int i = 0; i < n; ++i) {
        if (!visited[i]) {
            pq.emplace(_I->Node(i)->demand(), i);
        }
    }

    while (!pq.empty()) {
        int node = pq.top().second;
        double demand = pq.top().first;
        pq.pop();
        double best_increase = std::numeric_limits<double>::max();
        int best_route = -1;
        int best_pos = -1;

        // Try to insert into all routes at all possible positions
        for (int k = 0; k < K; ++k) {
            if (route_caps[k] < demand) continue;
            std::vector<int>& route = routes[k];
            // Try all possible insertion positions except depots
            for (size_t pos = 1; pos < route.size(); ++pos) {
                int prev = route[pos - 1];
                int next = route[pos];
                double old_cost = _I->cost(prev, next);
                double new_cost = _I->cost(prev, node) + _I->cost(node, next);
                double increase = new_cost - old_cost;
                if (increase < best_increase) {
                    best_increase = increase;
                    best_route = k;
                    best_pos = static_cast<int>(pos);
                }
            }
        }

        if (best_route != -1 && best_pos != -1) {
            // Insert node into the best position in the best route
            std::vector<int>& route = routes[best_route];
            int prev = route[best_pos - 1];
            int next = route[best_pos];
            route.insert(route.begin() + best_pos, node);
            route_caps[best_route] -= demand;
            total_cost += _I->cost(prev, node) + _I->cost(node, next) - _I->cost(prev, next);
            visited[node] = true;
            if (total_cost >= _bestUB) {
                // If the cost is already worse than the best known solution, return empty solution
                return new Solution(_I);
            }
        } else {
            return new Solution(_I); // No feasible insertion found
        }
    }

    // Create the solution from the routes
    Solution* sol = new Solution(_I);
    for (int k = 0; k < K; ++k) {
        const std::vector<int>& route = routes[k];
        if (route.size() <= 2) continue; // Skip empty or depot-only routes
        for (size_t i = 0; i < route.size() - 1; ++i) {
            int from = route[i];
            int to = route[i + 1];
            sol->addEdge(from, to, 1.0);
        }
    }

    return sol;
}

std::pair<std::vector<NumVar*>, std::vector<NumVar*>> BBSolver::_splitVars(Solution* sol, int i) const
{
    std::vector<FlowEdge*> edges;
    for (int j = 0; j < _I->nnodes(); j++) {
        if (j == i || _I->Node(j)->isInitialDepot()) continue; // Skip self-loops and depots
		double flow = sol->getFlow(i, j);
        FlowEdge* e = new FlowEdge(j, i, j, flow);
        edges.push_back(e);
        
    }

    sort(edges.begin(), edges.end(), [](FlowEdge* a, FlowEdge* b) {
        return a->getFlow() > b->getFlow();
        });


    // Check that edges has at least two elements and the first two have flow > 0
    if (edges.size() < 2 || edges[0]->getFlow() <= 0.0 || edges[1]->getFlow() <= 0.0) {
        // Clean up edges before returning
        for (FlowEdge* e : edges) {
            delete e;
        }
        throw std::runtime_error("Not enough valid edges with positive flow for splitting.");
    }

    // Split edges into left and right based on even/odd index
    std::vector<NumVar*> left;
    std::vector<NumVar*> right;
    for (size_t idx = 0; idx < edges.size(); ++idx) {
		NumVar* var = _LP->getVar(edges[idx]);
        if (idx % 2 == 0) {
			left.push_back(var);
        }
        else {
			right.push_back(var);
        }
    }

	// Clean up edges
    for (FlowEdge* e : edges) {
        delete e;
    }

	// Return the split variables
	return std::pair<std::vector<NumVar*>, std::vector<NumVar*>>(left, right);


}
