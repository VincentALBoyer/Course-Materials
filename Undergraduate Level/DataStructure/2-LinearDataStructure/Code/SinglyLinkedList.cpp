/*
 * SinglyLinkedList - A generic singly linked list implementation in C++.
 *
 * Template Parameters:
 *   T - The type of elements stored in the list.
 *
 * Features:
 *   - push_front(const T& value): Inserts an element at the front of the list.
 *   - push_back(const T& value): Inserts an element at the back of the list.
 *   - insert(Iterator pos, const T& value): Inserts an element after the node pointed by the iterator 'pos'.
 *   - pop_front(): Removes the element at the front of the list.
 *   - pop_back(): Removes the element at the back of the list.
 *   - empty() const: Checks if the list is empty.
 *   - size() const: Returns the number of elements in the list.
 *   - clear(): Removes all elements from the list.
 *   - front(): Returns a reference to the first element.
 *   - back(): Returns a reference to the last element.
 *   - print() const: Prints all elements in the list to std::cout.
 *   - begin(), end(): Returns iterators to the beginning and end of the list.
 *
 * Nested Types:
 *   - Iterator: A forward iterator for traversing the list.
 *     - operator*(): Dereferences to the element data.
 *     - operator++(): Advances the iterator.
 *     - operator==, operator!=: Compares iterators.
 *
 * Usage Example:
 *   SinglyLinkedList<int> list;
 *   list.push_back(1);
 *   list.push_front(0);
 *   for (auto it = list.begin(); it != list.end(); ++it) {
 *       std::cout << *it << " ";
 *   }
 *
 * Notes:
 *   - The Iterator is only valid for the list it was created from.
 *   - Insertion with insert() places the new element after the given iterator position.
 *   - The list manages its own memory and cleans up nodes in the destructor.
 */
#include <iostream>

template <typename T>
class SinglyLinkedList {
private:
    struct Node {
        T data;
        Node* next;
        Node(const T& value) : data(value), next(nullptr) {}
    };
    Node* head;
    size_t list_size;

public:
    class Iterator {
        Node* current;

        
        const SinglyLinkedList* parent; // Allow access to the list for validation

        // Make SinglyLinkedList a friend to allow access to private members
        friend class SinglyLinkedList;

        // Private constructor to ensure Iterator can only be created by SinglyLinkedList
        Iterator(Node* node, const SinglyLinkedList* list) : current(node), parent(list) {}

    public:

        T& operator*() const { return current->data; }
        Iterator& operator++() { current = current->next; return *this; }
        bool operator!=(const Iterator& other) const { return current != other.current || parent!=other.parent ; }
        bool operator==(const Iterator& other) const { return current == other.current && parent==other.parent ; }
    };

public:
    SinglyLinkedList() : head(nullptr), list_size(0) {}

    ~SinglyLinkedList() {
        clear();
    }

    void push_front(const T& value) {
        Node* new_node = new Node(value);
        new_node->next = head;
        head = new_node;
        ++list_size;
    }

    void push_back(const T& value) {
        Node* new_node = new Node(value);
        if (!head) {
            head = new_node;
        } else {
            Node* curr = head;
            while (curr->next) curr = curr->next;
            curr->next = new_node;
        }
        ++list_size;
    }

    Iterator insert(Iterator pos, const T& value) {
        if (pos == end()) {
            push_back(value);
            // Find the last node to return its iterator
            Node* curr = head;
            while (curr->next) curr = curr->next;
            return Iterator(curr,this);
        }
        Node* curr = pos.current;
        // Node* curr = head;
        // while (curr && Iterator(curr) != pos) {
        //     curr = curr->next;
        // }
        // if (!curr) {
        //     throw std::out_of_range("Iterator out of range");
        // }
        Node* new_node = new Node(value);
        new_node->next = curr->next;
        curr->next = new_node;
        ++list_size;
        return Iterator(new_node,this);
    }

    void pop_front() {
        if (!head) return;
        Node* temp = head;
        head = head->next;
        delete temp;
        --list_size;
    }

    void pop_back() {
        if (!head) return;
        if (!head->next) {
            delete head;
            head = nullptr;
        } else {
            Node* curr = head;
            while (curr->next->next) curr = curr->next;
            delete curr->next;
            curr->next = nullptr;
        }
        --list_size;
    }

    bool empty() const {
        return head == nullptr;
    }

    size_t size() const {
        return list_size;
    }

    void clear() {
        while (head) {
            pop_front();
        }
    }

    T& front() {
        if (!head) throw std::out_of_range("List is empty");
        return head->data;
    }

    T& back() {
        if (!head) throw std::out_of_range("List is empty");
        Node* curr = head;
        while (curr->next) curr = curr->next;
        return curr->data;
    }

    void print() const {
        Node* curr = head;
        while (curr) {
            std::cout << curr->data << " ";
            curr = curr->next;
        }
        std::cout << std::endl;
    }

    Iterator begin() { return Iterator(head, this);}
    Iterator end() { return Iterator(nullptr, this);}
};


int main() {
    SinglyLinkedList<int> list;

    list.push_back(1);
    list.push_back(2);
    list.push_front(0);
    list.print(); // Output: 0 1 2

    // Instead of exposing head, let's provide an iterator for safe traversal.
    std::cout << "Iterating over list: ";
    for (SinglyLinkedList<int>::Iterator it = list.begin(); it != list.end(); ++it) {
        std::cout << *it << " ";
    }
    std::cout << std::endl;

    list.pop_front();
    list.print(); // Output: 1 2

    list.pop_back();
    list.print(); // Output: 1

    list.insert(list.begin(), -1);
    list.print(); // Output: -1 1

    auto it = list.begin();
    ++it; // Points to 1
    list.insert(it, 0);
    list.print(); // Output: -1 0 1

    std::cout << "Front: " << list.front() << std::endl; // Output: Front: 1
    std::cout << "Back: " << list.back() << std::endl;   // Output: Back: 1

    return 0;
}