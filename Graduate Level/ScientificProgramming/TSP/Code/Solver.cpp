#include "Solver.hpp"

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
        output << setw(width) << "#_nodes"
            << setw(width) << " Y"
            <<setw(width)<<"Solution" 
           << endl;
        output << " ";
    }

    Solution* S = recoversolution();
    string fitness = "n/a";
    string feas = "n/a";
    if (S != NULL) {
        fitness = to_string(S->fitness());
        if (S->isFeasible()) feas = "yes";
        else feas = "no";
    }

    std::time_t end_time = std::chrono::system_clock::to_time_t(_end);
    std::tm now;
    localtime_s(&now, &end_time);

    string date = PadNumber(now.tm_mday) + "/" + PadNumber(now.tm_mon + 1) + "/" + to_string(now.tm_year + 1900);
    string time = PadNumber(now.tm_hour) + ":" + PadNumber(now.tm_min);

    output << setw(20) << date << setw(20) << time << setw(40) << RedString(_I->instancename(),40) << setw(width) << RedString(_name,width) << setw(width) << fitness << setw(width) << gap() << setw(width) << cpuTime() << setw(width) << feas ;
    output << setw(width) << _I->nnodes();

    int sumy = 0;
    if (S != NULL) {
        for (auto i : *S) {
            sumy++;
        }
    }
    output << setw(width) << sumy;

    //int counter = 0;
    if(S!=NULL){
        string sol="(";
        
        for(auto i : *S){
            sol+=to_string(i);
            if(i!=S->back()) sol+=",";
            //counter++;
        }
        sol+=")";
        
        output<<sol;
    }
    else output<<setw(width)<<"n/a";
    output << endl;
    
    delete S;

    output.close();
}