/*
    BTree.cpp

    This file implements a basic B-Tree data structure with insertion and search functionalities.
    The B-Tree is a balanced tree data structure commonly used in databases and filesystems to
    maintain sorted data and allow searches, sequential access, insertions, and deletions in logarithmic time.

    Classes:
    --------
    1. BTreeNode
        - Represents a node in the B-Tree.
        - Members:
            - keys: Vector of keys stored in the node.
            - children: Vector of pointers to child nodes.
            - leaf: Boolean indicating if the node is a leaf.
        - Methods:
            - findKey(int k): Finds the first key greater than or equal to k.
            - insertNonFull(int k): Inserts a key into a node that is not full.
            - splitChild(int i, BTreeNode* y): Splits the child node y at index i.
            - print(int level = 0): Prints the subtree rooted at this node in in-order traversal.
            - search(int k): Searches for a key in the subtree rooted at this node.

    2. BTree
        - Represents the B-Tree itself.
        - Members:
            - root: Pointer to the root node.
        - Methods:
            - insert(int k): Inserts a key into the B-Tree.
            - print(): Prints the entire B-Tree in in-order traversal.
            - search(int k): Searches for a key in the B-Tree.

    Usage:
    ------
    - The main function demonstrates insertion of several keys into the B-Tree,
      prints the tree structure, and performs search operations.

    Note:
    -----
    - The minimum degree (T) is set to 3, meaning each node (except root) must have at least T-1 keys
      and at most 2*T-1 keys.
    - This implementation does not include deletion logic.
*/
// BTree.cpp
// Basic B-Tree node structure (no full insertion logic)
#include <iostream>
#include <vector>
using namespace std;


class BTreeNode {
public:
    vector<int> keys;
    vector<BTreeNode*> children;
    bool leaf;
    unsigned int T; // Minimum degree of the B-Tree

    BTreeNode(bool leaf, unsigned int T) : leaf(leaf), T(T) {}

    // Find the first key greater than or equal to k
    int findKey(int k) {
        int idx = 0;
        while (idx < keys.size() && keys[idx] < k)
            ++idx;
        return idx;
    }

    // Insert a new key in a non-full node
    void insertNonFull(int k) {
        int i = keys.size() - 1;
        if (leaf) {
            keys.push_back(0);
            while (i >= 0 && keys[i] > k) {
                keys[i + 1] = keys[i];
                i--;
            }
            keys[i + 1] = k;
        } else {
            while (i >= 0 && keys[i] > k)
                i--;
            if (children[i + 1]->keys.size() == 2 * T - 1) {
                splitChild(i + 1, children[i + 1]);
                if (keys[i + 1] < k)
                    i++;
            }
            children[i + 1]->insertNonFull(k);
        }
    }

    // Split the child y of this node at index i
    void splitChild(int i, BTreeNode* y) {
        BTreeNode* z = new BTreeNode(y->leaf, T);
        for (int j = 0; j < T - 1; j++)
            z->keys.push_back(y->keys[j + T]);
        if (!y->leaf) {
            for (int j = 0; j < T; j++)
                z->children.push_back(y->children[j + T]);
        }
        y->keys.resize(T - 1);
        if (!y->leaf)
            y->children.resize(T);
        children.insert(children.begin() + i + 1, z);
        keys.insert(keys.begin() + i, y->keys[T - 1]);
    }

    // Print the tree in-order
    void print(int level = 0) {
        int i;
        for (i = 0; i < keys.size(); i++) {
            if (!leaf)
                children[i]->print(level + 1);
            for (int l = 0; l < level; ++l) cout << "  ";
            cout << keys[i] << endl;
        }
        if (!leaf)
            children[i]->print(level + 1);
    }

    // Search for a key in the subtree rooted with this node
    BTreeNode* search(int k) {
        int i = 0;
        while (i < keys.size() && k > keys[i])
            i++;
        if (i < keys.size() && keys[i] == k)
            return this;
        if (leaf)
            return nullptr;
        return children[i]->search(k);
    }
};

class BTree {
    BTreeNode* root;
    unsigned int T = 3; // Minimum degree (T) of the B-Tree
public:
    BTree(unsigned int t = 3) : T(t) { root = new BTreeNode(true, T); }

    void insert(int k) {
        if (root->keys.size() == 2 * T - 1) {
            BTreeNode* s = new BTreeNode(false, T);
            s->children.push_back(root);
            s->splitChild(0, root);
            int i = 0;
            if (s->keys[0] < k)
                i++;
            s->children[i]->insertNonFull(k);
            root = s;
        } else {
            root->insertNonFull(k);
        }
    }

    void print() {
        if (root)
            root->print();
    }

    bool search(int k) {
        return root->search(k) != nullptr;
    }
};

int main() {
    BTree tree;
    vector<int> values = {10, 20, 5, 6, 12, 30, 7, 17};
    for (int v : values)
        tree.insert(v);

    cout << "B-Tree structure (in-order):" << endl;
    tree.print();

    int searchKey = 6;
    cout << "Searching for " << searchKey << ": "
         << (tree.search(searchKey) ? "Found" : "Not found") << endl;

    searchKey = 15;
    cout << "Searching for " << searchKey << ": "
         << (tree.search(searchKey) ? "Found" : "Not found") << endl;

    return 0;
}
