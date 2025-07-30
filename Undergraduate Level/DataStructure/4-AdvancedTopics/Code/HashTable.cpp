// HashTable.cpp
// Simple hash table with chaining (integer key, string value)
#include <iostream>
#include <list>
#include <vector>
using namespace std;

class HashTable {
    static const int SIZE = 10;
    vector<list<pair<int, string>>> table;
public:
    HashTable() : table(SIZE) {}
    int hash(int key) { return key % SIZE; }
    void insert(int key, const string& value) {
        table[hash(key)].push_back({key, value});
    }
    bool search(int key, string& value) {
        for (auto& kv : table[hash(key)]) {
            if (kv.first == key) { value = kv.second; return true; }
        }
        return false;
    }
    void print() {
        for (int i = 0; i < SIZE; ++i) {
            cout << i << ": ";
            for (auto& kv : table[i]) cout << "(" << kv.first << ", " << kv.second << ") ";
            cout << endl;
        }
    }
};

int main() {
    HashTable ht;
    ht.insert(1, "apple");
    ht.insert(11, "banana");
    ht.insert(21, "cherry");
    ht.print();
    string val;
    if (ht.search(11, val)) cout << "Found: " << val << endl;
    else cout << "Not found" << endl;
    return 0;
}
