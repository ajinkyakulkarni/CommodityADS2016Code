#ifndef LCP_TREAP_H
#define LCP_TREAP_H
#include "../RandomNumberGeneration/Random.h"
#include "../Utils/GCFreeList.h"
#include "Treap.h"

namespace igmdk{

template<typename KEY, typename VALUE, typename INDEXED_COMPARATOR =
    LexicographicComparator<KEY> > class LcpTreap
{
    INDEXED_COMPARATOR comparator;
    struct Node
    {
        KEY key;
        VALUE value;
        Node *left, *right;
        unsigned int priority;
        unsigned short predLcp, succLcp;
        Node(KEY const& theKey, VALUE const& theValue): key(theKey),
            value(theValue), left(0), right(0), priority(GlobalRNG.next()),
            predLcp(0), succLcp(0) {}
    }* root;
    Freelist<Node> freelist;
    Node* rotateRight(Node* node)
    {
        Node* goingUp = node->left;
        node->predLcp = goingUp->succLcp;
        goingUp->succLcp = min(node->succLcp, goingUp->succLcp);
        node->left = goingUp->right;
        goingUp->right = node;
        return goingUp;
    }
    Node* rotateLeft(Node* node)
    {
        Node* goingUp = node->right;
        node->succLcp = goingUp->predLcp;
        goingUp->predLcp = min(node->predLcp, goingUp->predLcp);
        node->right = goingUp->left;
        goingUp->left = node;
        return goingUp;
    }
    int findLCP(KEY const& key, Node* node, int predM, int& m)
    {
        int lcp = predM == m ? node->predLcp : node->succLcp;
        if(m <= lcp)
        {
            while(m < comparator.getSize(key) &&
               comparator.isEqual(key, node->key, m)) ++m;
            lcp = m;
        }
        return lcp;
    }
    Node* constructFrom(Node* node)
    {
        Node* tree = 0;
        if(node)
        {
            tree = new(freelist.allocate())Node(node->key, node->value);
            tree->priority = node->priority;
            tree->predLcp = node->predLcp;
            tree->succLcp = node->succLcp;
            tree->left = constructFrom(node->left);
            tree->right = constructFrom(node->right);
        }
        return tree;
    }
public:
    typedef Node NodeType;
    LcpTreap(INDEXED_COMPARATOR theComparator = INDEXED_COMPARATOR()):
        root(0), comparator(theComparator) {}
    LcpTreap(LcpTreap const& other): comparator(other.comparator)
        {root = constructFrom(other.root);}
    LcpTreap& operator=(LcpTreap const&rhs){return genericAssign(*this,rhs);}

    Node** findPointer(KEY const& key)
    {
        Node** pointer = &root, *node;
        int m = 0, predM = 0;
        while(node = *pointer)
        {
            int lcp = findLCP(key, node, predM, m);
            if(comparator.getSize(key) == lcp &&
                comparator.getSize(node->key) == lcp) break;
            if(comparator.isLess(key, node->key, lcp)) pointer = &node->left;
            else
            {
                pointer = &node->right;
                predM = lcp;
            }
        }
        return pointer;
    }
    Node* findNode(KEY const& key){return *findPointer(key);}
    VALUE* find(KEY const& key)
    {
        Node* node = findNode(key);
        return node ? &node->value : 0;
    }
    Node* insertNode(Node* newNode, Node* node, int m)
    {
        if(!node) return newNode;
        int lcp = findLCP(newNode->key, node, newNode->predLcp, m);
        if(comparator.getSize(node->key) == lcp &&
           comparator.getSize(newNode->key) == lcp)
        {
            node->value = newNode->value;
            freelist.remove(newNode);
        }
        else
        {
            bool goLeft = comparator.isLess(newNode->key, node->key, lcp);
            (goLeft ? newNode->succLcp : newNode->predLcp) = lcp;
            Node*& chosenChild = goLeft ? node->left : node->right;
            chosenChild = insertNode(newNode, chosenChild, m);
            if(chosenChild->priority < node->priority)
                node = goLeft ? rotateRight(node) : rotateLeft(node);
        }
        return node;
    }
    void insert(KEY const& key, VALUE const& value)
    {
        root = insertNode(new(freelist.allocate())Node(key, value), root, 0);
    }
    Node* removeFound(Node* node)
    {
        Node *left = node->left, *right = node->right;
        if(left && right)
        {
            bool goRight = left->priority < right->priority;
            node = goRight ? rotateRight(node) : rotateLeft(node);
            Node*& child = goRight ? node->right : node->left;
            child = removeFound(child);
        }
        else
        {
            freelist.remove(node);
            node = left ? left : right;
        }
        return node;
    }
    void remove(KEY const& key)
    {
        Node** node = findPointer(key);
        if(node) *node = removeFound(*node);
    }
    NodeType* predecessor(KEY const& key)
    {
        int m = 0, predM = 0;
        for(Node* node = root; node;)
        {
            int lcp = findLCP(key, node, predM, m);
            if(!comparator.isLess(node->key, key, lcp)) node = node->left;
            else
            {
                if(!node->right) return node;
                node = node->right;
                predM = lcp;
            }
        }
        return 0;
    }
    NodeType* successor(KEY const& key)
    {
        int m = 0, predM = 0;
        for(Node* node = root; node;)
        {
            int lcp = findLCP(key, node, predM, m);
            if(!comparator.isLess(key, node->key, lcp))
            {
                node = node->right;
                predM = lcp;
            }
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
    NodeType* prefixSuccessor(KEY const& key)
    {
        int m = 0, predM = 0;
        for(Node* node = root; node;)
        {
            int lcp = findLCP(key, node, predM, m);
            if(comparator.getSize(key) == lcp ||
                !comparator.isLess(key, node->key, lcp))
            {
                node = node->right;
                predM = lcp;
            }
            else if(!node->left) return node;
            else node = node->left;
        }
        return 0;
    }
    typedef TreeIterator<Node> Iterator;
    Iterator end(){return Iterator(0);}
    Iterator buildIterator(Node* theNode)
    {
        if(!theNode) return end();
        Iterator i(theNode);
        int m = 0, predM = 0;
        for(Node* node = root; node && node != theNode;)
        {
            i.ancestors.push(node);
            int lcp = findLCP(theNode->key, node, predM, m);
            if(comparator.isLess(theNode->key, node->key, lcp))
                node = node->left;
            else
            {
                node = node->right;
                predM = lcp;
            }
        }
        return i;
    }
    Iterator begin(){return Iterator(root);}
};

}//end namespace
#endif
