#pragma once


#include <iostream>
#include <time.h> 
#include <iomanip>
#include <windows.h>
#include <fstream>
#include <string>

using namespace std;

enum State { ALIVE, DEAD };

struct Grid {
    int size;
    State** states;
};

Grid initEmptyGrid(int n);

void printGrid(const Grid& G);

void deleteGrid(Grid& G);

void randomizeStates(Grid& G, double p = 0.25);

void setPattern(Grid& G, string File);

int countAlives(const Grid& G, int i, int j);

void upgradeGridState(Grid& G);

Grid copyGrid(const Grid& G);