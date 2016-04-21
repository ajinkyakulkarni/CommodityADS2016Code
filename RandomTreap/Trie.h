#ifndef TRIE_H
#define TRIE_H
#include "../Utils/GCFreeList.h"
#include "../Utils/Vector.h"
#include "../RandomNumberGeneration/Random.h"
namespace igmdk{

template<typename ITEM, typename KEY_OBJECT = unsigned char, typename
    COMPARATOR = DefaultComparator<ITEM> > class TernaryTreapTrie
{
    COMPARATOR comparator;
    struct Node
    {
        KEY_OBJECT pivot;
        unsigned int priority;
        Node *next, *left, *right;
        ITEM* item;
        Node(KEY_OBJECT const& thePivot): next(0), left(0), right(0),
            item(0), pivot(thePivot), priority(GlobalRNG.next()){}
    }* root;
    Freelist<ITEM> itemFreelist;
    Freelist<Node> nodeFreelist;
    Node* rotateRight(Node* tree)
    {
        Node* goingUp = tree->left;
        tree->left = goingUp->right;
        goingUp->right = tree;
        return goingUp;
    }
    Node* rotateLeft(Node* tree)
    {
        Node* goingUp = tree->right;
        tree->right = goingUp->left;
        goingUp->left = tree;
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
    template<typename ACTION> void forEachNode(Node* node, ACTION& action)
    {
        if(node)
        {
            action(node);
            forEachNode(node->left);
            forEachNode(node->next);
            forEachNode(node->right);
        }
    }
    struct CopyAction
    {
        Vector<ITEM*>& result;
        CopyAction(Vector<ITEM*>& theResult): result(theResult){}
        void operator()(Node* node){if(node->item)result.append(node->item);}
    };
    Node* constructFrom(Node* node)
    {
        Node* result = 0;
        if(node)
        {
            result = new(nodeFreelist.allocate())Node(node->pivot);
            result->priority = node->priority;
            if(node->item)
                result->item = new(itemFreelist.allocate())ITEM(node->item);
            result->left = constructFrom(node->left);
            result->next = constructFrom(node->next);
            result->right = constructFrom(node->right);
        }
        return result;
    }
public:
    TernaryTreapTrie(COMPARATOR const& theComparator = COMPARATOR()):
        root(0), comparator(theComparator){}
    TernaryTreapTrie(TernaryTreapTrie const& other):
        comparator(other.comparator){root = constructFrom(other.root);}
    TernaryTreapTrie& operator=(TernaryTreapTrie const& rhs)
        {return genericAssign(*this, rhs);}

    Node* insertNode(Node* node, KEY_OBJECT* key, int keySize,
        ITEM const& item, int i)
    {
        if(!node) node = new(nodeFreelist.allocate())Node(key[i]);
        if(comparator.isEqual(key[i], node->pivot))
            if(i == keySize - 1)
            {
                if(node->item) *node->item = item;
                else node->item = new(itemFreelist.allocate())ITEM(item);
            }
            else node->next = insertNode(node->next, key, keySize, item,
                i + 1);
        else
        {
            bool goLeft = comparator.isLess(key[i], node->pivot);
            Node*& chosenChild = goLeft ? node->left : node->right;
            chosenChild = insertNode(chosenChild, key, keySize, item, i);
            if(chosenChild->priority < node->priority)
                node = goLeft ? rotateRight(node) : rotateLeft(node);
        }
        return node;
    }
    void insert(unsigned char* key, int keySize, ITEM const& item)
    {
        assert(keySize > 0);
        root = insertNode(root, key, keySize, item, 0);
    }
    Vector<ITEM*> prefixFind(KEY_OBJECT* key, int lcp)
    {
        Vector<ITEM*> result;
        CopyAction action(result);
        forEachNode(findNode(key, lcp), action);
        return result;
    }
    ITEM* longestMatch(KEY_OBJECT* key, int keySize)
    {
        Node* node = root, *result = 0;
        for(int i = 0; node && i < keySize;)
            if(comparator.isEqual(key[i], node->pivot))
            {
                result = node;
                if(i == keySize - 1) break;
                else {node = node->next; ++i;}
            }
            else if(comparator.isLess(key[i], node->pivot))
                node = node->left;
            else node = node->right;
        return result ? result->item : 0;
    }
    struct Handle
    {
        Node* node;
        int i;
        Handle():node(0),i(0){}
    };
    Node* findNodeIncremental(KEY_OBJECT* key, int keySize,Node*& node,
        int& i)
    {
        while(node && i < keySize)
        {
            if(comparator.isEqual(key[i], node->pivot))
            {
                if(i == keySize - 1) return node;
                else{node = node->next; ++i;}
            }
            else if(comparator.isLess(key[i], node->pivot))node = node->left;
            else node = node->right;
        }
        i = 0;
        node = 0;
        return 0;
    }
    ITEM* findIncremental(KEY_OBJECT* key, int keySize, Handle& h)
    {
        if(!h.node) h.node = root;
        Node* result = findNodeIncremental(key, keySize, h.node, h.i);
        return result ? result->item : 0;
    }
    Node* findNode(KEY_OBJECT* key, int keySize)
    {
        Node* temp = root;
        int i = 0;
        return findNodeIncremental(key, keySize, temp, i);
    }
    ITEM* find(KEY_OBJECT* key, int keySize)
    {
        Node* result = findNode(key, keySize);
        return result ? result->item : 0;
    }
    Node* removeR(Node* node, KEY_OBJECT* key, int keySize, int i)
    {
        if(node)
        {
            bool isEqual = comparator.isEqual(key[i], node->pivot);
            if(isEqual && i == keySize - 1)
            {//remove found item
                if(!node->item) return node;
                itemFreelist.remove(node->item);
                node->item = 0;
            }
            else
            {//go to next node
                Node** child;
                if(isEqual){child = &node->next; ++i;}
                else if(comparator.isLess(key[i], node->pivot))
                    child = &node->left;
                else child = &node->right;
                *child = removeR(*child, key, keySize, i);
            }
            if(!node->item && !node->next)
            {//remove empty node
                Node* left = node->left, *right = node->right;
                nodeFreelist.remove(node);
                node = (left || right) ? join(left, right) : 0;
            }
        }
        return node;
    }
    void remove(KEY_OBJECT* key, int keySize)
    {
        assert(keySize > 0);
        root = removeR(root, key, keySize, 0);
    }
};
}
#endif
