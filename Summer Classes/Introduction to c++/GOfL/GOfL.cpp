// GOfL.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <time.h> 
#include <iomanip>

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

int main()
{
    int size = 10;

    Grid G = initEmptyGrid(size);

    randomizeStates(G, 0.8);

    printGrid(G);

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
    int w = 5;
    for (int i = 0; i < G.size; i++) {
        for (int j = 0; j < G.size; j++) {
            if (G.states[i][j] == State::ALIVE)
                cout << setw(2) << (char)254u;
            else cout << setw(2) << '_';
        }

        cout << endl;
    }
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


// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
