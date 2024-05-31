// Author: Mario Talevski
#include <iostream>
#include <iomanip>
#include <Windows.h>

using namespace std;

constexpr int BLACK_COLOR = 0;
constexpr int BLUE_COLOR = 9;
constexpr int GREEN_COLOR = 10;
constexpr int RED_COLOR = 12;
constexpr int WHITE_COLOR = 15;

enum State{ALIVE, DEAD};

struct Grid {
	int gen;
	int size;
	State** states;
};

Grid initEmptyGrid(int size);

void randomizeStates(Grid& G, double p = 0.1);

Grid copyGrid(const Grid& G);

void deleteGrid(Grid& G);

void resetCursorPosition();

void printGrid(const Grid& G);

void updateGridState(Grid& G);

int countAlives(const Grid& G, int i, int j);

int main(int argc, char argv[]) {

	Grid G = initEmptyGrid(50);

	randomizeStates(G, 0.5);

	//G.states[1][1] = State::ALIVE;
	//G.states[1][3] = State::ALIVE;
	//G.states[2][2] = State::ALIVE;
	//G.states[2][3] = State::ALIVE;
	//G.states[3][2] = State::ALIVE;


	while (G.gen++ < 500) {
		printGrid(G);
		updateGridState(G);
		Sleep(100);
	}
	

	return 0;
}

void printGrid(const Grid& G)
{
	resetCursorPosition();
	//system("CLS");
	int w = 2;
	HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
	for (int i = 0; i < G.size; ++i) {
		for (int j = 0; j < G.size; ++j) {
			if (G.states[i][j] == State::DEAD) 
				SetConsoleTextAttribute(handle, BLACK_COLOR);
			else 
				SetConsoleTextAttribute(handle, BLUE_COLOR);
			
			cout << setw(w) << (char)254u;
		}
		cout << endl;
	}
	SetConsoleTextAttribute(handle, WHITE_COLOR);

	//Counting the number of alive cells
	int n = 0;
	for (int i = 0; i < G.size; ++i)
		for (int j = 0; j < G.size; ++j) n += (G.states[i][j] == State::ALIVE);


	int W = 10;
	cout << setw(W) << "Gen: " << setw(W) << G.gen << endl;
	cout << setw(W) << "Alive: " << setw(W) << n << endl;

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

void randomizeStates(Grid& G, double p)
{
	srand(time(0));
	for (int i = 0; i < G.size; ++i) {
		for (int j = 0; j < G.size; ++j) {
			double r = rand() % 1001 / 1000.0;
			if (r < p) G.states[i][j] = State::ALIVE;
			else G.states[i][j] = State::DEAD;
		}
	}
}

Grid copyGrid(const Grid& G)
{
	Grid H = initEmptyGrid(G.size);
	for (int i = 0; i < G.size; ++i) {
		for (int j = 0; j < G.size; ++j) H.states[i][j] = G.states[i][j];
	}
	return H;
}

void deleteGrid(Grid& G)
{
	for (int i = 0; i < G.size; ++i) delete[] G.states[i];
	delete[] G.states;
}

void resetCursorPosition() {
	COORD topLeft = { 0, 0 };
	HANDLE console = GetStdHandle(STD_OUTPUT_HANDLE);
	//CONSOLE_SCREEN_BUFFER_INFO screen;
	//DWORD written;

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
