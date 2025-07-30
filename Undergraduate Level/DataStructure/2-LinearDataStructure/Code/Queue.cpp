/**
 * This file defines a template class `Queue` that provides a queue (FIFO) interface
 * for any data type. The queue is implemented using a singly linked list, with
 * pointers to the front and rear nodes, and supports standard queue operations.
 *
 * @tparam T Type of elements stored in the queue.
 *
 * @class Queue
 * @brief A generic queue supporting enqueue, dequeue, peek, and size operations.
 *
 * Public Methods:
 * - Queue() : Constructs an empty queue.
 * - ~Queue() : Destructor. Deallocates all nodes in the queue.
 * - void enqueue(const T& value) : Adds an element to the rear of the queue.
 * - void dequeue() : Removes the element at the front of the queue. Throws std::out_of_range if the queue is empty.
 * - T& peek() : Returns a reference to the front element. Throws std::out_of_range if the queue is empty.
 * - bool isEmpty() const : Checks if the queue is empty.
 * - size_t size() const : Returns the number of elements in the queue.
 *
 * Example usage:
 * @code
 * Queue<int> q;
 * q.enqueue(10);
 * q.enqueue(20);
 * q.enqueue(30);
 * while (!q.isEmpty()) {
 *     std::cout << q.peek() << " ";
 *     q.dequeue();
 * }
 * @endcode
 */
#include <iostream>

template <typename T>
class Queue {
private:
    struct Node {
        T data;
        Node* next;
        Node(const T& val) : data(val), next(nullptr) {}
    };
    Node* front;
    Node* rear;
    size_t count;

public:
    Queue() : front(nullptr), rear(nullptr), count(0) {}

    ~Queue() {
        while (!isEmpty()) {
            dequeue();
        }
    }

    void enqueue(const T& value) {
        Node* newNode = new Node(value);
        if (isEmpty()) {
            front = rear = newNode;
        } else {
            rear->next = newNode;
            rear = newNode;
        }
        ++count;
    }

    void dequeue() {
        if (isEmpty()) {
            throw std::out_of_range("Queue is empty");
        }
        Node* temp = front;
        front = front->next;
        delete temp;
        --count;
        if (front == nullptr) {
            rear = nullptr;
        }
    }

    T& peek() {
        if (isEmpty()) {
            throw std::out_of_range("Queue is empty");
        }
        return front->data;
    }

    bool isEmpty() const {
        return front == nullptr;
    }

    size_t size() const {
        return count;
    }
};

// Example usage

int main() {
    Queue<int> q;
    q.enqueue(10);
    q.enqueue(20);
    q.enqueue(30);
    while (!q.isEmpty()) {
        std::cout << q.peek() << " ";
        q.dequeue();
    }
    return 0;
}
