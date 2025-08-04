#pragma once
#include "ORUtils.hpp"


enum class DistanceMode {Euclidian=0, Pseudoeuclidian=1, Fischetti=2, Geographical=3, Explicit=4};

enum class NodeType {InitialDepot=0, FinalDepot=1, Node=3};

class Instance;

class node{
    int _id;
    double _x;
    double _y;
	double _demand;
    NodeType _type;

    friend class Instance;
    
public:

    node(int i=0, double x=0, double y=0, double demand=0, NodeType t=NodeType::Node): _id(i), _x(x), _y(y), _demand(demand), _type(t) {};
    ~node() {};
    
    int id() const {return _id;};
    double x() const {return _x;};
    double y() const {return _y;};
	double demand() const { return _demand; };  
    NodeType type() const {return _type;};
    bool isDepot() {return _type == NodeType::InitialDepot || _type == NodeType::FinalDepot;};
    bool isInitialDepot() const {return _type == NodeType::InitialDepot;};
    bool isFinalDepot() const {return _type == NodeType::FinalDepot;};
    bool isNode() const {return _type == NodeType::Node;};

    void print();


    
};

class Instance {
    
    string _name;

    string _knownsolpath;
    
    DistanceMode _DM;

    double** _D;        //Distance matrix for explicit instances
    
    int _nnodes;
    
    int _InitialDepotId;
    
    int _FinalDepotId;

	double _Q; 		//Vehicle capacity  

	int _p;      //Number of vehicles
    
    node** _Nodes;



    void readExplicit(ifstream& input, string edge_weight_format);

    void setTravellingTime(int i, int j, double t) {
        if(i==j) return;    //No self loops
        int ii=max(i,j);
        int jj=min(i,j);
        _D[ii][jj] = t;
    }
    


public:

    bool collinear(const node* a, const node* b, const node* c) {
        // Implement the collinear function
        return (b->_y - a->_y) * (c->_x - b->_x) == (c->_y - b->_y) * (b->_x - a->_x);
    }

	Instance(string filename);
	~Instance();
    
	string name() const { return _name; };

    string knownsolpath() const { return _knownsolpath; }

    int indexmap(int realid) {
		if (realid > 0 && _Nodes[realid - 1]->id() == realid) {
			return realid - 1;
		}
        else{
			auto it = find_if(_Nodes, _Nodes + _nnodes, [realid](node* n) {return n->id() == realid; });
			if (it == _Nodes + _nnodes) {
				throw runtime_error("Node " + to_string(realid) + " not found");
			}
			return (int)distance(_Nodes, it);
		}
	}

    int nnodes() const {return _nnodes;}
    
    node* InitialDepot() {return _Nodes[_InitialDepotId];};
    node* FinalDepot() {return _Nodes[_FinalDepotId];};
	int InitialDepotId() const { return _InitialDepotId; };
	int FinalDepotId() const { return _FinalDepotId; };
    node* Node(int i) {return _Nodes[i];}
    double VehicleCap() const { return _Q; }
    int FleetSize() const { return _p; }
    
    double travellingtime(node* Ni, node* Nj);
	double travellingtime(int i, int j) { return travellingtime(Node(i), Node(j)); }
	double cost(int i, int j) { return travellingtime(i, j); }

	double length(vector<int>& tour);

    double demand(const std::vector<int>& S) const {
        return accumulate(S.begin(), S.end(), 0.0, [&](double sum, int i) {
            return sum + _Nodes[i]->demand();
            });
    }
        
	void print();
 
};