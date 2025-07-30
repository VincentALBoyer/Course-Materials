#include <iostream>
#include <vector>

class TreeNode {
public:
    int value;
    std::vector<TreeNode*> children;

    TreeNode(int val) : value(val) {}
    ~TreeNode() {
        for (auto child : children) {
            delete child;
        }
    }
};

class SimpleTree {
private:
    TreeNode* root;

    void print(TreeNode* node, int depth) const {
        if (!node) return;
        // Print current level nodes
        for (int i = 0; i < depth; ++i) std::cout << "  ";
        std::cout << node->value << " -> {";
        for (auto child : node->children) {
            std::cout << " " << child->value;
        }
        std::cout << " }"<< std::endl;
        // Recurse for each child
        for (auto child : node->children) {
            print(child, depth + 1);
        }
    }

    TreeNode* find(TreeNode* node, int val) const {
        if (!node) return nullptr;
        if (node->value == val) return node;
        for (auto child : node->children) {
            TreeNode* res = find(child, val);
            if (res) return res;
        }
        return nullptr;
    }

public:
    SimpleTree(int rootVal) {
        root = new TreeNode(rootVal);
    }

    ~SimpleTree() {
        delete root;
    }

    // Add child to parent node with value parentVal
    bool addChild(int parentVal, int childVal) {
        TreeNode* parent = find(root, parentVal);
        if (!parent) return false;
        parent->children.push_back(new TreeNode(childVal));
        return true;
    }

    // Print tree
    void print() const {
        print(root, 0);
    }

    // Find if a value exists in the tree
    bool contains(int val) const {
        return find(root, val) != nullptr;
    }
};


int main() {
    SimpleTree tree(1);
    tree.addChild(1, 2);
    tree.addChild(1, 3);
    tree.addChild(2, 4);
    tree.addChild(2, 5);
    tree.addChild(3, 6);

    std::cout << "Tree structure:" << std::endl;
    tree.print();

    int searchVal = 5;
    std::cout << "Contains " << searchVal << "? " 
              << (tree.contains(searchVal) ? "Yes" : "No") << std::endl;

    searchVal = 7;
    std::cout << "Contains " << searchVal << "? " 
              << (tree.contains(searchVal) ? "Yes" : "No") << std::endl;

    return 0;
}