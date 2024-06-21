// GOfL.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "GameOfLife.h"


int main(int argc, const char* argv[])
{

    cout << "argc: " << argc << endl;
    for (int i = 0; i < argc; i++)
        cout << string(argv[i]) << endl;

    string File = "";
    if(argc>=2) File=argv[1];

    int size = 50;
    if(argc>=3) size = atoi(argv[2]);

    int ngen = 500;
    if (argc >= 4) ngen = atoi(argv[3]);

    Grid G = initEmptyGrid(size);

    if (File == "")
        randomizeStates(G, 0.5);
    else
        setPattern(G, File);

    printGrid(G);

    for (int i = 0; i < ngen; i++) {
        //Print Grid
        printGrid(G);

        //Update Grid State
        upgradeGridState(G);

        //Wait 
        Sleep(100);
    }
    

    deleteGrid(G);




}




// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
