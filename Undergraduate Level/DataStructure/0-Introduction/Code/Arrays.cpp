/*
    This program demonstrates dynamic memory allocation for 1D and 2D arrays in C++.

    Features:
    - Dynamically allocates a 1D array of integers, initializes it with multiples of 10, prints its contents, and deallocates the memory.
    - Dynamically allocates a 2D array of integers (as an array of pointers), initializes it with sequential values, prints its contents in matrix form, and properly deallocates all allocated memory.

    Key Concepts:
    - Usage of `new` and `delete[]` for dynamic memory management.
    - Proper initialization and cleanup of dynamically allocated arrays to prevent memory leaks.
    - Nested loops for handling 2D array initialization and output.
*/

#include <iostream>

int main() {
    // Dynamic allocation of 1D array
    int n = 5;
    int* arr1D = new int[n];
    for (int i = 0; i < n; ++i) {
        arr1D[i] = i * 10;
    }
    std::cout << "1D Array: ";
    for (int i = 0; i < n; ++i) {
        std::cout << arr1D[i] << " ";
    }
    std::cout << std::endl;
    delete[] arr1D;

    // Dynamic allocation of 2D array
    int rows = 3, cols = 4;
    int** arr2D = new int*[rows];
    for (int i = 0; i < rows; ++i) {
        arr2D[i] = new int[cols];
    }
    // Initialize and print 2D array
    std::cout << "2D Array:" << std::endl;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            arr2D[i][j] = i * cols + j;
            std::cout << arr2D[i][j] << " ";
        }
        std::cout << std::endl;
    }
    // Deallocate 2D array
    for (int i = 0; i < rows; ++i) {
        delete[] arr2D[i];
    }
    delete[] arr2D;

    return 0;
}