#include "Solution.hpp"



Solution::Solution(string filename, Instance* I) : vector<int>(), _I(I) {
   
    ifstream file(filename);

    if(!file.is_open()) {
        throw invalid_argument("Error: could not open file " + filename);
    }

    //Check if the file contains a NAME and DIMENSION
    string name = Instance::getData(file, "NAME");
    string dimension = Instance::getData(file, "DIMENSION");

    if(name == "n/a") {
        throw runtime_error("Error: could not find keyword NAME in file " + filename);
    }
    if(dimension == "n/a") {
        throw runtime_error("Error: could not find keyword DIMENSION in file " + filename);
    }
    cout << "Name: " << name << endl;
    cout << "Dimension: " << dimension << endl;

    // Check if the instance and the solution match
    int n=stoi(dimension);
    if(n != I->nnodes() || name.find(I->instancename()) == string::npos) {
        throw runtime_error("Error: instance and solution do not match");
    }

    // Check if the file contains a TOUR_SECTION
    if(!Instance::gotoSection(file, "TOUR_SECTION")) {
        throw invalid_argument("Error: could not find section TOUR_SECTION in file " + filename);
    }

    // Read the solution
    for (int i = 0; i < _I->nnodes(); i++) {
        int node;
        file >> node;
        push_back(node-1);
    }
    file.close();
}