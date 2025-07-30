/**
 * This file defines a Stack class that supports basic stack operations such as push, pop, and peek.
 * The stack dynamically resizes its underlying array when it becomes full.
 *
 * Classes:
 *  - Stack: Provides methods to manipulate a stack of integers.
 *
 * Stack Public Methods:
 *  - Stack(int size_ = INITIAL_SIZE): Constructor to initialize the stack with a given or default size.
 *  - ~Stack(): Destructor to free allocated memory.
 *  - bool isEmpty() const: Checks if the stack is empty.
 *  - bool isFull() const: Checks if the stack is full.
 *  - void push(int value): Pushes an integer onto the stack, resizing if necessary.
 *  - int pop(): Removes and returns the top element of the stack. Prints a message and returns -1 if empty.
 *  - int peek() const: Returns the top element without removing it. Prints a message and returns -1 if empty.
 *
 * Example Usage:
 *  Stack stack;
 *  stack.push(10);
 *  int top = stack.peek();
 *  int popped = stack.pop();
 */
#include <iostream>


class Stack {
private:
    static const int INITIAL_SIZE = 100;

    int* arr;
    int top;
    int capacity;
public:
    Stack(int size_ = INITIAL_SIZE) : top(-1), capacity(size_) {
        arr = new int[capacity];
    }

    ~Stack() {
        delete[] arr;
    }

    bool isEmpty() const {
        return top == -1;
    }

    bool isFull() const {
        return top == capacity - 1;
    }

    void push(int value) {
        if (isFull()) {
            // Increase stack size (dynamic resizing)
            int newCapacity = capacity * 2;
            int* newArr = new int[newCapacity];
            for (int i = 0; i < capacity; ++i) {
                newArr[i] = arr[i];
            }
            delete[] arr;
            arr = newArr;
            capacity = newCapacity;
        }
        arr[++top] = value;
    }

    int pop() {
        if (isEmpty()) {
            std::cout << "Stack Underflow\n";
            return -1;
        }
        return arr[top--];
    }

    int peek() const {
        if (isEmpty()) {
            std::cout << "Stack is empty\n";
            return -1;
        }
        return arr[top];
    }
};

int main() {
    Stack stack;

    stack.push(10);
    stack.push(20);
    stack.push(30);

    std::cout << "Top element is: " << stack.peek() << std::endl;

    std::cout << "Popped element: " << stack.pop() << std::endl;
    std::cout << "Top element after pop: " << stack.peek() << std::endl;

    // Popping all elements
    while (!stack.isEmpty()) {
        std::cout << "Popped element: " << stack.pop() << std::endl;
    }

    return 0;
}
