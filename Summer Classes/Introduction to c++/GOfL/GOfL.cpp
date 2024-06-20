// GOfL.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <time.h> 
#include <iomanip>
#include <windows.h>
#include <fstream>
#include <string>

using namespace std;

enum State{ALIVE, DEAD};

struct Grid{
    int size;
    State** states;
};

Grid initEmptyGrid(int n);

void printGrid(const Grid& G);

void deleteGrid(Grid& G);

void randomizeStates(Grid& G, double p=0.25);

void setPattern(Grid& G, string File);

int countAlives(const Grid& G, int i, int j);

void upgradeGridState(Grid& G);

Grid copyGrid(const Grid& G);

int main()
{
    int size = 50;
    int ngen = 500;

    Grid G = initEmptyGrid(size);

    //randomizeStates(G, 0.5);

    setPattern(G, "D:/GitHub/Course-Materials/Summer Classes/Introduction to c++/Patterns/fun.pgol");

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

Grid initEmptyGrid(int n)
{
    Grid G;
    G.size = n;
    G.states = new State*[n];
    for (int i = 0; i < n; i++) {
        G.states[i] = new State[n];
        for (int j = 0; j < n; j++)
            G.states[i][j] = State::DEAD;
    }

    return G;
}

void printGrid(const Grid& G)
{
    COORD topLeft = { 0,0 };
    HANDLE console = GetStdHandle(STD_OUTPUT_HANDLE);
    SetConsoleCursorPosition(console, topLeft);
    
    SetConsoleTextAttribute(console, 12);
    int w = 5;
    for (int i = 0; i < G.size; i++) {
        for (int j = 0; j < G.size; j++) {
            if (G.states[i][j] == State::ALIVE)
                cout << setw(2) << (char)254u;
            else cout << setw(2) << ' ';
        }

        cout << endl;
    }
    SetConsoleTextAttribute(console, 15);
}

void deleteGrid(Grid& G)
{
    for (int i = 0; i < G.size; i++) delete[] G.states[i];
    delete[] G.states;
}

void randomizeStates(Grid& G, double p)
{
    srand(time(NULL));
    for (int i = 0; i < G.size; i++) {
        for (int j = 0; j < G.size; j++) {
            double r = (double)(rand() % 1001) / 1000.0;
            if (r < p)
                G.states[i][j] = State::ALIVE;
            else G.states[i][j] = State::DEAD;
        }
    }
}

void setPattern(Grid& G, string File)
{
    ifstream input(File.c_str());

    if (!input.is_open()) {
        cout << "File " << File << " not found!" << endl;
        exit(25);
    }

    streampos pos = input.tellg();

    int col_dim = 0; //for index j
    while (input.get() != '\n') col_dim++;

    input.seekg(pos);

    int row_dim = 0;
    string line;
    while (getline(input, line)) 
        row_dim++;

    input.clear();
    input.seekg(pos);

    int x = (int)floor(G.size / 2.0 - col_dim/2.0);
    int y = (int)floor(G.size / 2.0 - row_dim/2.0);

    for (int i = 0; i < row_dim; i++) {
        for (int j = 0; j < col_dim+1; j++) {
            char c = input.get();
            if (c == '*' || c == 'X')
                G.states[x + i][y + j] = ALIVE;
            //cout << c;
        }
        //cout << endl;
    }
    //cout << endl;
    

       

}

int countAlives(const Grid& G, int i, int j)
{
    int nalives = 0;

    for (int k = -1; k <= 1; k++) {
        for (int l = -1; l <= 1; l++) if(k!=0 || l!=0){
            int x = (k + i + G.size) % G.size;
            int y = (l + j + G.size) % G.size;
            if (G.states[x][y] == ALIVE) nalives++;
        }
    }

    return nalives;
}

void upgradeGridState(Grid& G)
{
    Grid tempG = copyGrid(G);
    for (int i = 0; i < G.size; i++)  for (int j = 0; j < G.size; j++) {
        int n = countAlives(tempG, i, j);
        if (tempG.states[i][j] == ALIVE && (n == 2 || n == 3)) 
            G.states[i][j] = ALIVE;
        else if (tempG.states[i][j] == DEAD && n == 3) 
            G.states[i][j] = ALIVE;
        else 
            G.states[i][j] = DEAD;
    }
    
    deleteGrid(tempG);
}

Grid copyGrid(const Grid& G)
{
    Grid H = initEmptyGrid(G.size);
    for (int i = 0; i < G.size; i++) {
        for (int j = 0; j < G.size; j++) {
            H.states[i][j] = G.states[i][j];
        }
    }


    return H;
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
