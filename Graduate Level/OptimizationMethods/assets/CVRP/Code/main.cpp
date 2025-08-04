#include "Instance.hpp"
#include "LPSolvers/LPSolver.hpp"
#include "IPSolvers/IPSolver.hpp"






#include <filesystem>
namespace fs = std::filesystem;


using namespace std;

// Test routine to solve all .vrp files in a folder
void test_solve_folder(const string& folder, int timlim, Solver* (*solver_factory)(Instance*)) {
    vector<string> vrp_files;
    for (const auto& entry : fs::directory_iterator(folder)) {
        if (entry.is_regular_file() && entry.path().extension() == ".vrp") {
            vrp_files.push_back(entry.path().string());
        }
    }

    // Select up to 10 random .vrp files from the folder
    std::shuffle(vrp_files.begin(), vrp_files.end(), ORRandom::engine());
    if (vrp_files.size() > 10) {
        vrp_files.resize(10);
    }
    for (const auto& filename : vrp_files) {
        cout << "Solving: " << filename << endl;
        Instance* I = nullptr;
        Solver* solver = nullptr;
        try {
            I = new Instance(filename);
            solver = solver_factory(I);
            solver->setparam(Solver::TimLim, timlim);
            solver->solve();
           
			solver->save("test_log.txt");
        } catch (const exception& e) {
            cout << "Error solving " << filename << ": " << e.what() << endl;
			exit(1);
        }
        delete solver;
        delete I;

		
    }
    exit(100); // Exit after solving one instance for testing purposes
}

int main(int argc, const char * argv[]) {

    //test_solve_folder("D:\\GitHub\\ResearchProjects\\CVRP\\Instances\\SetA\\", 1800, [](Instance* I) -> Solver* { return new CuttingPlaneSolver(I); });

    string filename;
    int timlim=1800;
    if (argc>1) {
        filename=argv[1];
        if (argc>2) {
            timlim=atoi(argv[2]);
        }
    } else {
        cout << "Please provide the filename of the instance to solve" << endl;
        return 0;
    }

    Instance* I = NULL;
    try{
        I=new Instance(filename);
        I->print();
    } catch (exception& e) {
        cout << "Error: " << e.what() << endl;
        return 1;
    }

    Env env(I);
    //exit(1000);

    LPSolver* LP = new CuttingPlaneSolver(I);
    LP->setparam(Solver::TimLim, timlim);
    try {
        LP->solve();
    }
    catch (const std::exception& e) {
        std::cout << "Error during LP->solve(): " << e.what() << std::endl;
        delete LP;
        delete I;
        ORRandom::destroy();
        return 1;
    }
    Solver* IP = new BBSolver(LP);
 //   LPSolver* LP = new CuttingPlaneSolver(I); 
	//LP->setparam(Solver::TimLim, timlim);
 //   try {
    //    LP->solve();
    //} catch (const std::exception& e) {
    //    std::cout << "Error during LP->solve(): " << e.what() << std::endl;
    //    delete LP;
    //    delete I;
    //    ORRandom::destroy();
    //    return 1;
    //}
 //   exit(100);
	//Solver* IP = new BBSolver(LP);

    IP->setparam(Solver::TimLim, timlim - LP->cpuTime());

    try {
        IP->solve();
    } catch (const std::exception& e) {
        std::cout << "Error during IP->solve(): " << e.what() << std::endl;
        delete IP;
        delete LP;
        delete I;
        ORRandom::destroy();
        return 1;
    }
    Solution* sol = IP->recoversolution();

    if (sol != NULL) {
        sol->print();
        sol->isFeasible();
    }
	else {
		cout << "No solution found" << endl;
	}

    try{    
        string tourfile = ORUtils::extractpath(filename) + "\\" + ORUtils::extractfilename(filename) + ".sol";
        Solution* sol2=new Solution(tourfile,I);
		cout << "Known solution: " << sol2->fitness() << endl;
        sol2->print();

        if (sol != NULL) {
            cout << "Feasible solution: " << sol->fitness() << endl;
            cout << "Gap: " << 100 * (sol->fitness() - sol2->fitness()) / sol2->fitness() << "%" << endl;
        }
        cout << "Optimal solution: " << sol2->fitness() << endl;
        

        delete sol2;
    }
    catch (exception& e) {
        cout << e.what() << endl;
    }

	//cout << "Gap: " << S->gap() << endl;
	//cout << "Time: " << S->cpuTime() << endl;
    
    
    if(sol!=NULL) delete sol;
	delete IP;
	delete LP;
    delete I;

	ORRandom::destroy();

    return 0;

}