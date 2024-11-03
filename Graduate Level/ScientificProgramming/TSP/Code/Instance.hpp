#pragma once
#include "Header.hpp"


enum class DistanceMode {Euclidian=0, Pseudoeuclidian=1, Fischetti=2, Geographical=3, Explicit=4};

enum class NodeType {InitialDepot=0, FinalDepot=1, Node=3};

class Instance;

class node{
    int _id;
    double _x;
    double _y;
    NodeType _type;

    friend class Instance;
    
public:

    node(int i=0, double x=0, double y=0, NodeType t=NodeType::Node): _id(i), _x(x), _y(y), _type(t) {};
    ~node() {};
    
    int id() {return _id;};
    double x() const {return _x;};
    double y() const {return _y;};
    NodeType type() {return _type;};
    bool isDepot() {return _type == NodeType::InitialDepot || _type == NodeType::FinalDepot;};
    bool isInitialDepot() {return _type == NodeType::InitialDepot;};
    bool isFinalDepot() {return _type == NodeType::FinalDepot;};
    bool isNode() {return _type == NodeType::Node;};

    void print();
    
};

class Instance {
    
    string _name;
    
    DistanceMode _DM;

    double** _D;        //Distance matrix for explicit instances
    
    int _nnodes;
    
    int _InitialDepotId;
    
    int _FinalDepotId;
    
    node** _Nodes;

    static double coordtogeo(double val) {
        double pi = 3.141592;
        double deg = (int)val;
        double min = val - deg;
        return pi * (deg + 5.0 * min / 3.0) / 180;
    }

    void readExplicit(ifstream& input, string edge_weight_format);

    void setTravellingTime(int i, int j, double t) {
        if(i==j) return;    //No self loops
        int ii=max(i,j);
        int jj=min(i,j);
        _D[ii][jj] = t;
    }
    
    // Trim from the start (left)
    static inline std::string ltrim(std::string s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
            return !std::isspace(ch);
        }));
        return s;
    }
    
    // Trim from the end (right)
    static inline std::string rtrim(std::string s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
            return !std::isspace(ch);
        }).base(), s.end());
        return s;
    }
    
    // Trim from both ends
    static inline std::string trim(std::string s) {
        return ltrim(rtrim(s));
    }

public:

    string static getData(ifstream& input, string keyword) {
        input.clear();
        input.seekg(0, ios::beg);
        string line;
        while (getline(input, line)) {
            size_t pos = line.find(keyword);
            if (pos != string::npos) {
                pos = line.find(':',pos);
                if (pos == string::npos) {
                    return "n/a";
                }
                string s=line.substr(pos + 1);
                return trim(s);
            }
        }
        return "n/a";
    }

    bool static gotoSection(ifstream& input, string keyword) {
        input.clear();
        input.seekg(0, ios::beg);
        string line;
        while (getline(input, line)) {
            if (line.find(keyword) != string::npos) {
                return true;
            }
        }
        return false;
    }


	Instance(string filename);
	~Instance();
    
    string instancename() {
        size_t pos = _name.find_last_of("/");
        return _name.substr(pos + 1);
        
    }

    int nnodes() {return _nnodes;}
    
    node* InitialDepot() {return _Nodes[_InitialDepotId];};
    node* FinalDepot() {return _Nodes[_FinalDepotId];};
    node* Node(int i) {return _Nodes[i];}
    
    double travellingtime(node* Ni, node* Nj);
        
	void print();
 
};