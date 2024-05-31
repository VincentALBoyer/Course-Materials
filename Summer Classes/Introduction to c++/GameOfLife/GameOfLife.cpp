// Author: Mario Talevski
#include <iostream>
#include <iomanip>
#include <Windows.h>

using namespace std;

enum State{ALIVE, DEAD};

struct Grid {
	int gen;
	int size;
	State** states;
};

Grid initEmptyGrid(int size);

void printGrid(const Grid& G);

void updateGridState(Grid& G);

int countAlives(const Grid& G, int i, int j);


int main(int argc, char argv[]) {

	Grid G = initEmptyGrid(20);

	printGrid(G);

	return 0;
}

void printGrid(const Grid& G)
{
	system("CLS");
	int w = 2;
	HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
	for (int i = 0; i < G.size; ++i) {
		for (int j = 0; j < G.size; ++j) {
			if (G.states[i][j] == State::DEAD) 
				SetConsoleTextAttribute(handle, 15);
			else 
				SetConsoleTextAttribute(handle, 10);
			
			cout << setw(w) << (char)254u;
		}
		cout << endl;
	}
	SetConsoleTextAttribute(handle, 15);
	cout << "Gen: " << G.gen << endl;

}

Grid initEmptyGrid(int size)
{
	Grid G;
	G.size = size;
	G.gen = 0;
	G.states = new State * [size];
	for (int i = 0; i < G.size; ++i) {
		G.states[i] = new State[G.size];
		for (int j = 0; j < G.size; ++j) G.states[i][j] = State::DEAD;
	}
	return G;
}
