// Heaps.cpp
// Simple min heap implementation
#include <iostream>
#include <vector>
using namespace std;

class MinHeap {
    vector<int> heap;
    void heapifyUp(int i) {
        while (i > 0 && heap[(i-1)/2] > heap[i]) {
            swap(heap[i], heap[(i-1)/2]);
            i = (i-1)/2;
        }
    }
    void heapifyDown(int i) {
        int n = heap.size();
        while (2*i+1 < n) {
            int j = 2*i+1;
            if (j+1 < n && heap[j+1] < heap[j]) j++;
            if (heap[i] <= heap[j]) break;
            swap(heap[i], heap[j]);
            i = j;
        }
    }
public:
    void insert(int val) {
        heap.push_back(val);
        heapifyUp(heap.size()-1);
    }
    int extractMin() {
        int minVal = heap[0];
        heap[0] = heap.back();
        heap.pop_back();
        heapifyDown(0);
        return minVal;
    }
    void print() {
        for (int v : heap) cout << v << " ";
        cout << endl;
    }
};

int main() {
    MinHeap h;
    h.insert(5); h.insert(3); h.insert(8); h.insert(1);
    h.print();
    cout << "Extracted min: " << h.extractMin() << endl;
    h.print();
    return 0;
}
