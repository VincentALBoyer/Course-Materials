/**
 * @file Functions.cpp
 * @brief Demonstrates basic function usage in C++ including declaration, definition, and calling.
 *
 * This file contains examples of:
 * - Declaring and defining functions in C++
 * - Calling functions from the main function
 * - Performing simple operations such as greeting, addition, and squaring a number
 *
 * Functions included:
 * - void greet(): Prints a welcome message to the console.
 * - int add(int a, int b): Returns the sum of two integers.
 * - double square(double num): Returns the square of a given double value.
 *
 */

#include <iostream>

// Function declaration
void greet(); 
int add(int a, int b); 
double square(double num); 

int main() {
    // Calling functions
    greet();
    
    int sum = add(5, 3);
    std::cout << "Sum of 5 and 3: " << sum << std::endl;

    double result = square(4.5);
    std::cout << "Square of 4.5: " << result << std::endl;

    return 0;
}

// Function definition for greeting
void greet() {
    std::cout << "Hello! Welcome to C++ Functions.\n";
}

// Function definition for addition
int add(int a, int b) {
    return a + b;
}

// Function definition for squaring a number
double square(double num) {
    return num * num;
}