/*
 * Demonstrates basic usage of pointers in C++.
 *
 * This program:
 * - Declares an integer variable `num` and initializes it.
 * - Declares a pointer `ptr` that stores the address of `num`.
 * - Prints the value of `num`, its address, the value stored in `ptr`, and the value pointed to by `ptr`.
 * - Modifies the value of `num` through the pointer and displays the updated value.
 *
 * Key concepts:
 * - Variable declaration and initialization
 * - Pointer declaration and assignment
 * - Dereferencing pointers to access or modify the value at a memory address
 * - Outputting variable values and addresses using `std::cout`
 */

#include <iostream>

int main() {
    // Declare an integer variable
    int num = 42;

    // Declare a pointer to store the address of num
    int* ptr = &num;

    // Display variable, address, and pointer value
    std::cout << "Value of num: " << num << std::endl;
    std::cout << "Address of num: " << &num << std::endl;
    std::cout << "Value stored in ptr (address of num): " << ptr << std::endl;
    std::cout << "Value pointed to by ptr: " << *ptr << std::endl;

    // Modify the value using the pointer
    *ptr = 100;
    std::cout << "New value of num after modification through pointer: " << num << std::endl;
     
    return 0;
}