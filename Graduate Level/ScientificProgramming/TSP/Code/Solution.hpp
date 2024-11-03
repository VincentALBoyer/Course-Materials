#include "Instance.hpp"


class Solution : public vector<int> {
    
    Instance* _I;

public:
    Solution(Instance* I): vector<int>() , _I(I) {};
    Solution(string filename, Instance* I);
    ~Solution() {};
    
    double fitness() {
        double f = 0;
        for (int i = 0; i < _I->nnodes(); i++) {
            f += _I->travellingtime(_I->Node((*this)[i]), _I->Node((*this)[(i + 1) % _I->nnodes()]));
        }
        return f;
    }
    
    bool isFeasible() {
        vector<string> Error;
        if (size() != _I->nnodes()) {
            Error.push_back("Number of nodes is different from the number of nodes in the instance");
        }
        for (int i = 0; i < _I->nnodes(); i++) {
            if (find(begin(), end(), i) == end()) {
                Error.push_back("Node " + to_string(i) + " is missing");
            }
        }
        if (Error.size() > 0) {
            cout << "Solution is not feasible" << endl;
            for (auto e : Error) {
                cout << e << endl;
            }
            return false;
        }
        else return true;
    }
    
    void print() {
        cout << "Solution: ";
        for (auto i : *this) {
            cout << _I->Node(i)->id() << " ";
        }
        cout << endl;
    }

    bool containsUndirectedEdge(int i, int j) {
        for (int k = 0; k < _I->nnodes(); k++) {
            if ((*this)[k] == i && (*this)[(k + 1) % _I->nnodes()] == j) return true;
            if ((*this)[k] == j && (*this)[(k + 1) % _I->nnodes()] == i) return true;
        }
        return false;
    }
    
};