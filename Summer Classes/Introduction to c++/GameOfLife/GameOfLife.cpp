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

Grid copyGrid(const Grid& G);

void deleteGrid(Grid& G);

void printGrid(const Grid& G);

void updateGridState(Grid& G);

int countAlives(const Grid& G, int i, int j);

void clear() {
	COORD topLeft = { 0, 0 };
	HANDLE console = GetStdHandle(STD_OUTPUT_HANDLE);
	CONSOLE_SCREEN_BUFFER_INFO screen;
	DWORD written;

	//GetConsoleScreenBufferInfo(console, &screen);
	//FillConsoleOutputCharacterA(
	//	console, ' ', screen.dwSize.X * screen.dwSize.Y, topLeft, &written
	//);
	//FillConsoleOutputAttribute(
	//	console, FOREGROUND_GREEN | FOREGROUND_RED | FOREGROUND_BLUE,
	//	screen.dwSize.X * screen.dwSize.Y, topLeft, &written
	//);
	SetConsoleCursorPosition(console, topLeft);
}


int main(int argc, char argv[]) {

	Grid G = initEmptyGrid(20);

	G.states[1][1] = State::ALIVE;
	G.states[1][3] = State::ALIVE;
	G.states[2][2] = State::ALIVE;
	G.states[2][3] = State::ALIVE;
	G.states[3][2] = State::ALIVE;


	while (G.gen++ < 100) {
		printGrid(G);
		updateGridState(G);
		Sleep(100);
	}
	

	return 0;
}

void printGrid(const Grid& G)
{
	clear();
	//system("CLS");
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

void updateGridState(Grid& G)
{
	Grid tempG = copyGrid(G);
	for (int i = 0; i < G.size; ++i) for (int j = 0; j < G.size; ++j){
		int n = countAlives(tempG, i, j);
		if (tempG.states[i][j] == State::ALIVE && (n == 2 || n == 3))
			G.states[i][j] = State::ALIVE;
		else if (tempG.states[i][j] == State::DEAD && n == 3)
			G.states[i][j] = State::ALIVE;
		else
			G.states[i][j] = State::DEAD;
	}

	deleteGrid(tempG);
}

int countAlives(const Grid& G, int i, int j)
{
	int nalives = 0;
	for (int k = -1; k <= 1; ++k) for (int l = -1; l <= 1; ++l) if(k!=0 || l!=0){
		
		int x = (k + i + G.size) % G.size;
		int y = (l + j + G.size) % G.size;
		nalives += (int)(G.states[x][y]==State::ALIVE);
	}
	return nalives;
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

Grid copyGrid(const Grid& G)
{
	Grid H = initEmptyGrid(G.size);
	for (int i = 0; i < G.size; ++i) {
		for (int j = 0; j < G.size; ++j) H.states[i][j]= G.states[i][j];
	}
	return H;
}

void deleteGrid(Grid& G)
{
	for (int i = 0; i < G.size; ++i) delete[] G.states[i];
	delete[] G.states;
}
