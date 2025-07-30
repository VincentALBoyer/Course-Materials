/*
 * Demonstrates static and dynamic array allocation in C++.
 *
 * - Allocates and initializes a static array of integers.
 * - Prints the elements of the static array.
 * - Dynamically allocates an array of integers using `new`.
 * - Initializes the dynamic array with values.
 * - Prints the elements of the dynamic array.
 * - Properly deallocates the dynamic array using `delete[]` to prevent memory leaks.
 */

#include <iostream>

int main() {
    // Static array allocation
    int staticArr[5] = {1, 2, 3, 4, 5};
    std::cout << "Static array elements: ";
    for (int i = 0; i < 5; ++i) {
        std::cout << staticArr[i] << " ";
    }
    std::cout << std::endl;

    // Dynamic array allocation
    int size = 5;
    int* dynamicArr = new int[size];
    for (int i = 0; i < size; ++i) {
        dynamicArr[i] = (i + 1) * 10;
    }
    std::cout << "Dynamic array elements: ";
    for (int i = 0; i < size; ++i) {
        std::cout << dynamicArr[i] << " ";
    }
    std::cout << std::endl;

    // Don't forget to free dynamically allocated memory
    delete[] dynamicArr;

    return 0;
}