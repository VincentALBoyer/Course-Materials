/*
    This program demonstrates the declaration, initialization, and usage of basic variable types in C++.

    Variable types included:
    - int:         Stores integer values (e.g., age)
    - float:       Stores floating-point numbers with single precision (e.g., pi)
    - double:      Stores floating-point numbers with double precision (e.g., precisePi)
    - char:        Stores single characters (e.g., grade)
    - bool:        Stores boolean values (true/false) (e.g., isStudent)
    - std::string: Stores sequences of characters (e.g., name)

    The program initializes each variable with a sample value and prints their values to the console.
*/

#include <iostream>
#include <string>  // For using strings
#include <limits>  // For std::numeric_limits

int main() {
    // Integer variable
    int age = 25;

    // Floating point variable
    float pi = 3.14;
    double precisePi = 3.1415926535; // Higher precision

    // Character variable
    char grade = 'A';

    // Boolean variable
    bool isStudent = true;

    // String variable
    std::string name = "Vincent";

    // Displaying the values
    std::cout << "Integer (age): " << age << std::endl;
    std::cout << "Float (pi): " << pi << std::endl;
    std::cout << "Double (precisePi): " << precisePi << std::endl;
    std::cout << "Character (grade): " << grade << std::endl;
    std::cout << "Boolean (isStudent): " << std::boolalpha << isStudent << std::endl;
    std::cout << "String (name): " << name << std::endl;

    // Variables Range
    std::cout << "\nVariable Ranges:" << std::endl;
    std::cout << "int range: " << std::numeric_limits<int>::min() << " to " << std::numeric_limits<int>::max() << std::endl;
    std::cout << "float range: " << std::numeric_limits<float>::min() << " to " << std::numeric_limits<float>::max() << std::endl;
    std::cout << "double range: " << std::numeric_limits<double>::min() << " to " << std::numeric_limits<double>::max() << std::endl;
    std::cout << "char range: " << static_cast<int>(std::numeric_limits<char>::min()) << " to " 
              << static_cast<int>(std::numeric_limits<char>::max()) << std::endl;
    std::cout << "bool range: " << std::boolalpha << std::numeric_limits<bool>::min() << " to "
              << std::numeric_limits<bool>::max() << std::endl;


    return 0;
}
