/*
 * LoopsAndConditions.cpp
 *
 * Demonstrates basic control flow constructs in C++:
 *  - if-else conditions: Checks if a number is positive, negative, or zero.
 *  - switch statement: Matches a number against specific cases.
 *  - for loop: Iterates a fixed number of times.
 *  - while loop: Performs countdown using a condition.
 *  - do-while loop: Executes a block at least once and repeats based on a condition.
 *
 * Output is printed to the console to illustrate the behavior of each construct.
 */

#include <iostream>

int main() {
    int num = 5;

    // If-else conditions
    std::cout << "Checking conditions:\n";
    if (num > 0) {
        std::cout << "num is positive.\n";
    } else if (num < 0) {
        std::cout << "num is negative.\n";
    } else {
        std::cout << "num is zero.\n";
    }

    // Switch statement example
    std::cout << "Checking switch case:\n";
    switch (num) {
        case 1:
            std::cout << "num is 1.\n";
            break;
        case 5:
            std::cout << "num is 5.\n";
            break;
        default:
            std::cout << "num is neither 1 nor 5.\n";
    }

    // For loop example
    std::cout << "Using a for loop:\n";
    for (int i = 1; i <= 5; i++) {
        std::cout << "Iteration " << i << "\n";
    }

    // While loop example
    std::cout << "Using a while loop:\n";
    int count = 5;
    while (count > 0) {
        std::cout << "Countdown: " << count << "\n";
        count--;
    }

    // Do-while loop example
    std::cout << "Using a do-while loop:\n";
    int x = 1;
    do {
        std::cout << "Value: " << x << "\n";
        x++;
    } while (x <= 3);

    return 0;
}