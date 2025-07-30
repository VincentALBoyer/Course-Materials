/*
    Inheritance.cpp

    This file demonstrates object-oriented programming concepts in C++,
    specifically inheritance, virtual functions, and multiple inheritance.

    Classes:
    --------
    1. Shape (Base class)
        - Represents a generic shape.
        - Declares a virtual draw() method for polymorphic behavior.
        - Declares a virtual destructor.

    2. Circle (Derived class, virtually inherits Shape)
        - Represents a circle.
        - Overrides draw() to provide circle-specific drawing behavior.

    3. Square (Derived class, virtually inherits Shape)
        - Represents a square.
        - Overrides draw() to provide square-specific drawing behavior.

    4. ColoredShape (Derived class, inherits from both Circle and Square)
        - Demonstrates multiple inheritance.
        - Provides drawBoth() to call draw() from both Circle and Square.
        - Overrides draw() to resolve ambiguity and provide custom behavior.

    Main Function:
    --------------
    - Instantiates a ColoredShape object.
    - Calls drawBoth() to demonstrate multiple inheritance.
    - Demonstrates polymorphism by calling draw() via a Shape pointer.

    Notes:
    ------
    - Virtual inheritance is used to avoid the diamond problem.
    - Overriding draw() in ColoredShape resolves ambiguity from multiple inheritance.
*/
#include <iostream>

// Base class
class Shape {
public:
    virtual void draw() const {
        std::cout << "Drawing a generic shape." << std::endl;
    }
    virtual ~Shape() = default;
};

// Derived class: Circle
class Circle : virtual public Shape {
public:
    void draw() const override {
        std::cout << "Drawing a circle." << std::endl;
    }
};

// Derived class: Square
class Square : virtual public Shape {
public:
    void draw() const override {
        std::cout << "Drawing a square." << std::endl;
    }
};

// Multiple inheritance: ColoredShape
class ColoredShape : public Circle, public Square {
public:
    // Calls draw() of both Circle and Square
    void drawBoth() const {
        std::cout << "[ColoredShape] ";
        Circle::draw();
        std::cout << "[ColoredShape] ";
        Square::draw();
    }
    // Optionally, override draw() to clarify ambiguity
    void draw() const override {
        std::cout << "Drawing a colored shape (circle and square):" << std::endl;
        drawBoth();
    }
};

int main() {
    ColoredShape cs;
    cs.drawBoth();

    // Optionally, demonstrate polymorphism
    Shape* s = &cs;
    s->draw(); // Calls ColoredShape::draw()

    return 0;
}