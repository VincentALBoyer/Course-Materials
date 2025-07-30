/*
    AbstractClass.cpp

    This file demonstrates the use of abstract classes and pure virtual functions in C++.

    Classes:
    --------
    1. Shape (Abstract Base Class)
        - Represents a generic geometric shape.
        - Contains two pure virtual functions:
            - double area() const: Computes and returns the area of the shape.
            - void draw() const: Outputs a textual representation of the shape.
        - Declares a virtual destructor for proper cleanup of derived objects.

    2. Circle (Derived from Shape)
        - Represents a circle with a given radius.
        - Implements:
            - double area() const: Calculates the area using πr².
            - void draw() const: Outputs a message describing the circle.

    3. Rectangle (Derived from Shape)
        - Represents a rectangle with a given width and height.
        - Implements:
            - double area() const: Calculates the area as width × height.
            - void draw() const: Outputs a message describing the rectangle.

    main():
    -------
    - Demonstrates polymorphism by creating Shape pointers to Circle and Rectangle objects.
    - Calls draw() and area() on each object to show dynamic binding.
    - Cleans up allocated memory.

    Usage:
    ------
    Compile and run to see how abstract classes and virtual functions enable polymorphic behavior in C++.
*/
#include <iostream>
using namespace std;

// Abstract base class
class Shape {
public:
    // Pure virtual function
    virtual double area() const = 0;
    virtual void draw() const = 0;
    virtual ~Shape() {}
};

class Circle : public Shape {
    double radius;
public:
    Circle(double r) : radius(r) {}
    double area() const override {
        return 3.14159 * radius * radius;
    }
    void draw() const override {
        cout << "Drawing a circle with radius " << radius << endl;
    }
};

class Rectangle : public Shape {
    double width, height;
public:
    Rectangle(double w, double h) : width(w), height(h) {}
    double area() const override {
        return width * height;
    }
    void draw() const override {
        cout << "Drawing a rectangle " << width << " x " << height << endl;
    }
};

int main() {
    Shape* s1 = new Circle(5.0);
    Shape* s2 = new Rectangle(4.0, 6.0);

    s1->draw();
    cout << "Area: " << s1->area() << endl;

    s2->draw();
    cout << "Area: " << s2->area() << endl;

    delete s1;
    delete s2;
    return 0;
}