#ifndef TREAP_H
#define TREAP_H
#include "../RandomNumberGeneration/Random.h"
#include "../Utils/GCFreeList.h"
#include "../Utils/Stack.h"
namespace igmdk{

template<typename NODE> class TreeIterator
{
    NODE* current;
public:
    Stack<NODE*> ancestors;
    TreeIterator(NODE* root){current = root;}
    TreeIterator& operator++()
    {
        assert(current);
        if(current->right)
        {
            ancestors.push(current);
            current = current->right;
            for(; current; current = current->left) ancestors.push(current);
            current = ancestors.pop();
        }
        else if(!ancestors.isEmpty())
        {
            NODE* parent = ancestors.pop();
            while(parent->right == current && !ancestors.isEmpty())
            {
                current = parent;
                parent = ancestors.pop();
            }
            current = parent->right == current ? 0 : parent;
        }
        else current = 0;
        return *this;
    }
    NODE& operator*(){assert(current); return *current;}
    NODE* operator->(){assert(current); return current;}
    bool operator!=(TreeIterator const& rhs){return current != rhs.current;}
};

template<typename KEY, typename VALUE, typename COMPARATOR =
    DefaultComparator<KEY> > class Treap
{
    COMPARATOR comparator;
    struct Node
    {
        KEY key;
        VALUE value;
        Node *left, *right;
        unsigned int priority;
        Node(KEY const& theKey, VALUE const& theValue): key(theKey),
        value(theValue), left(0), right(0), priority(GlobalRNG.next()){}
    }* root;
    Freelist<Node> freelist;
    Node* rotateRight(Node* node)
    {
        Node* goingUp = node->left;
        node->left = goingUp->right;
        goingUp->right = node;
        return goingUp;
    }
    Node* rotateLeft(Node* node)
    {
        Node* goingUp = node->right;
        node->right = goingUp->left;
        goingUp->left = node;
        return goingUp;
    }
    Node* join(Node* left, Node* right)
    {
        if(!left) return right;
        if(!right) return left;
        if(left->priority < right->priority)
        {
            left->right = join(left->right, right);
            return left;
        }
        else
        {
            right->left = join(left, right->left);
            return right;
        }
    }
    Node* constructFrom(Node* node)
    {
        Node* tree = 0;
        if(node)
        {
            tree = new(freelist.allocate())Node(node->key, node->value);
            tree->priority = node->priority;
            tree->left = constructFrom(node->left);
            tree->right = constructFrom(node->right);
        }
        return tree;
    }
public:
    typedef Node NodeType;
    bool isEmpty(){return !root;}
    Treap(COMPARATOR const& theComparator = COMPARATOR()):
        root(0), comparator(theComparator){}
    Treap(Treap const& other): comparator(other.comparator)
        {root = constructFrom(other.root);}
    Treap& operator=(Treap const& rhs){return genericAssign(*this, rhs);}

    Node** findPointer(KEY const& key)
    {
        Node *node, **pointer = &root;
        while((node = *pointer) && !comparator.isEqual(key, node->key))
            pointer = &(comparator.isLess(key, node->key) ?
                node->left : node->right);
        return pointer;
    }
    NodeType* findNode(KEY const& key){return *findPointer(key);}
    VALUE* find(KEY const& key)
    {
        Node* node = findNode(key);
        return node ? &node->value : 0;
    }

    NodeType* predecessor(KEY const& key)
    {
        for(Node* node = root; node;)
        {
            if(!comparator.isLess(node->key, key)) node = node->left;
            else if(!node->right) return node;
            else node = node->right;
        }
        return 0;
    }
    NodeType* successor(KEY const& key)
    {
        for(Node* node = root; node;)
        {
            if(!comparator.isLess(key, node->key)) node = node->right;
            else if(!node->left) return node;
            else node = node->left;
        }
        return 0;
    }
    NodeType* inclusiveSuccessor(KEY const& key)
    {
        Node* node = findNode(key);
        return node ? node : successor(key);
    }

    Node* insertInto(KEY const& key, VALUE const& value, Node* node)
    {
        if(!node) return new(freelist.allocate())Node(key, value);
        if(comparator.isEqual(key, node->key)) node->value = value;
        else
        {
            bool goLeft = comparator.isLess(key, node->key);
            Node*& chosenChild = goLeft ? node->left : node->right;
            chosenChild = insertInto(key, value, chosenChild);
            if(chosenChild->priority < node->priority)
                node = goLeft ? rotateRight(node) : rotateLeft(node);
        }
        return node;
    }
    void insert(KEY const& key, VALUE const& value)
        {root = insertInto(key, value, root);}
    void remove(KEY const& key)
    {
        Node **pointer = findPointer(key), *node = *pointer;
        if(node)
        {
            *pointer = join(node->left, node->right);
            freelist.remove(node);
        }
    }
    NodeType* findMin()
    {
        Node* node = root;
        if(node) while(node->left) node = node->left;
        return node;
    }
    NodeType* findMax()
    {
        Node* node = root;
        if(node) while(node->right) node = node->right;
        return node;
    }
    typedef TreeIterator<Node> Iterator;
    Iterator begin(){return buildIterator(findMin());}
    Iterator end(){return Iterator(0);}
    Iterator buildIterator(Node* theNode)
    {
        if(!theNode) return end();
        Iterator i(theNode);
        Node* node = root;
        while(node && !comparator.isEqual(theNode->key, node->key))
        {
            i.ancestors.push(node);
            node = comparator.isLess(theNode->key, node->key) ?
                node->left : node->right;
        }
        return i;
    }
};

}//end namespace
#endif
