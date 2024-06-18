// Array.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <time.h> 
#include <iomanip>

using namespace std;

bool** initGrid(int n);

void printGrid(bool** G, int n);

void deleteGrid(bool** G, int n);

int main()
{
    int size = 10;
   
    bool** Grid = initGrid(size);

    printGrid(Grid, size);

    deleteGrid(Grid, size);

  
   

}

bool** initGrid(int n)
{
    bool** G = new bool* [n];
    for (int i = 0; i < n; i++) {
        G[i] = new bool[n];
        for (int j = 0; j < n; j++)
            G[i][j] = false;
    }

    return G;
}

void printGrid(bool** G, int n)
{
    int w = 5;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (G[i][j] == false)
                cout << setw(2) << (char)254u;
        }

        cout << endl;
    }
}

void deleteGrid(bool** G, int n)
{
    for (int i = 0; i < n; i++) delete[] G[i];
    delete[] G;
}
