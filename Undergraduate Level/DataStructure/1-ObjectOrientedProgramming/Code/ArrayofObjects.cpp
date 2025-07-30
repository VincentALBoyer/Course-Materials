/*
    ArrayofObjects.cpp

    This program demonstrates the use of classes, dynamic arrays, and copy constructors in C++.

    Features:
    - Defines a `Student` class without an explicit copy constructor.
    - Defines a `Teacher` class with an explicit copy constructor that outputs a message when called.
    - Dynamically allocates arrays of `Student` and `Teacher` objects.
    - Assigns values to array elements and displays their contents.
    - Demonstrates the invocation of the copy constructor for the `Teacher` class.
    - Properly deallocates dynamically allocated memory.

    Key Concepts:
    - Default vs. explicit copy constructors.
    - Dynamic memory allocation for arrays of objects.
    - Object assignment and copy semantics in C++.
*/
#include <iostream>
#include <string>
using namespace std;

// Class without explicit copy constructor
class Student {
    string name;
    int age;
public:
    

    Student(string n = "", int a = 0) : name(n), age(a) {}

    void display() const {
        cout << "Name: " << name << ", Age: " << age << endl;
    }
};

// Class with explicit copy constructor
class Teacher {
    string name;
    int id;
public:
    

    Teacher(string n = "", int i = 0) : name(n), id(i) {}

    // Explicit copy constructor
    Teacher(const Teacher& t) {
        name = t.name;
        id = t.id;
        cout << "Copy constructor called for Teacher: " << name << endl;
    }

    void display() const {
        cout << "Teacher: " << name << ", ID: " << id << endl;
    }
};

int main() {
    int nStudents = 2;
    int nTeachers = 2;

    // Dynamic allocation for students (with default copy constructor for simple class)
    Student* students = new Student[nStudents];
    students[0] = Student("Alice", 20);
    students[1] = Student("Bob", 21);

    cout << "Students:" << endl;
    for (int i = 0; i < nStudents; ++i)
        students[i].display();

    // Dynamic allocation for teachers (with explicit copy constructor)
    Teacher* teachers = new Teacher[nTeachers];
    teachers[0] = Teacher("Mr. Smith", 101);
    teachers[1] = Teacher("Ms. Lee", 102);

    cout << "\nTeachers:" << endl;
    for (int i = 0; i < nTeachers; ++i)
        teachers[i].display();

    // Demonstrate copy constructor
    Teacher t2 = teachers[0]; // Copy constructor called
    t2.display();

    // Clean up
    delete[] students;
    delete[] teachers;

    return 0;
}