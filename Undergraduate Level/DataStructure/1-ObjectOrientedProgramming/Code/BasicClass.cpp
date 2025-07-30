
/*
    BasicClass.cpp

    This file demonstrates a simple C++ class named 'Student' that encapsulates
    basic object-oriented programming concepts such as constructors, destructors,
    private member variables, and public member functions.

    Classes:
    --------
    Student:
        - Represents a student with a name and an ID.
        - Constructor initializes the student's name and ID, and prints a creation message.
        - Destructor prints a destruction message when a Student object is destroyed.
        - display(): Outputs the student's name and ID to the console.

    Usage:
    ------
    The main function creates a Student object and displays its information.

    Example Output:
    ---------------
    Student object created: Alice
    Name: Alice, ID: 101
    Student object destroyed: Alice
*/
#include <iostream>
#include <string>

class Student {
private:
    std::string name;
    int id;
public:
    // Constructor
    Student(const std::string& studentName, int studentId)
        : name(studentName), id(studentId) {
        std::cout << "Student object created: " << name << std::endl;
    }

    // Destructor
    ~Student() {
        std::cout << "Student object destroyed: " << name << std::endl;
    }

    // Method to display student info
    void display() const {
        std::cout << "Name: " << name << ", ID: " << id << std::endl;
    }
};

int main() {
    Student s1("Alice", 101);
    s1.display();

    Student s2("Bob", 102);
    s2.display();

    return 0;
}