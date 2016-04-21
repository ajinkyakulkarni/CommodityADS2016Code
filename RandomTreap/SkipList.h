#ifndef SKIPLIST_H
#define SKIPLIST_H
#include "../RandomNumberGeneration/Random.h"
#include "../Utils/GCFreeList.h"
namespace igmdk{

template<typename KEY, typename VALUE, typename COMPARATOR =
    DefaultComparator<VALUE> > class SkipList
{
    COMPARATOR comparator;
    enum{MAX_HEIGHT = 32};
    struct Node
    {
        KEY key;
        VALUE value;
        Node** tower;
        Node(KEY const& theKey, VALUE const& theValue, int height):
            key(theKey), value(theValue), tower(new Node*[height]) {}
        ~Node(){delete[] tower;}
    }* head[MAX_HEIGHT];
    Freelist<Node> freelist;
    int currentLevel;
public:
    typedef Node NodeType;
    SkipList(COMPARATOR const& theComparator = COMPARATOR()):
        currentLevel(0), comparator(theComparator)
        {for(int i = 0; i < MAX_HEIGHT; ++i) head[i] = 0;}
    SkipList(SkipList const& rhs): currentLevel(0),
        comparator(rhs.comparator)
    {//order of items with nonunique keys in copy is reversed
        for(int i = 0; i < MAX_HEIGHT; ++i) head[i] = 0;
        for(Node* node = rhs.head[0]; node; node = node->tower[0])
            insert(node->key, node->value, false);
    }
    SkipList& operator=(SkipList const&rhs){return genericAssign(*this,rhs);}

    NodeType* predecessor(KEY const& key)
    {
        Node **tower = head, *pred = 0;
        for(int level = currentLevel; level >= 0; --level)
            for(Node* node; (node = tower[level]) && comparator.isLess(
                node->key, key); tower = node->tower) pred = node;
        return pred;
    }
    NodeType* inclusiveSuccessor(KEY const& key)
    {
        Node* pred = predecessor(key);
        return pred? pred->tower[0] : findMin();
    }
    NodeType* findNode(KEY const& key)
    {
        Node* node = inclusiveSuccessor(key);
        return node && !comparator.isLess(node->key, key) ? node : 0;
    }
    NodeType* successor(KEY const& key)
    {//is inclusiveSuccessor with non-unique keys
        Node* node = inclusiveSuccessor(key);
        if(node && !comparator.isLess(node->key, key)) node = node->tower[0];
        return node;
    }
    VALUE* find(KEY const& key)
    {
        Node* result = findNode(key);
        return result ? &result->value : 0;
    }

    Node* insert(KEY const& key, VALUE const& value, bool unique = true)
    {
        if(unique)
        {
            Node* result = findNode(key);
            if(result)
            {
                result->value = value;
                return result;
            }
        }
        int newLevel = min(MAX_HEIGHT - 1, GlobalRNG.geometric(0.632));
        Node* newNode = new(freelist.allocate())Node(key, value,
            newLevel + 1);
        if(currentLevel < newLevel) currentLevel = newLevel;
        Node** tower = head;
        for(int level = currentLevel; level >= 0; --level)
        {
            for(Node* node; (node = tower[level]) &&
                comparator.isLess(node->key, key); tower = node->tower);
            if(level <= newLevel)
            {
                newNode->tower[level] = tower[level];
                tower[level] = newNode;
            }
        }
        return newNode;
    }
    void remove(KEY const& key)
    {
        Node **tower = head, *result = 0;
        for(int level = currentLevel; level >= 0; --level)
            for(Node* node; (node = tower[level]) &&
                !comparator.isLess(key, node->key); tower = node->tower)
                if((result && result == node) ||
                    (!result && !comparator.isLess(node->key, key)))
                {//can simplify the "if" if don't allow nonunique keys
                    tower[level] = node->tower[level];
                    if(!head[currentLevel]) --currentLevel;
                    result = node;
                    break;
                }
        if(result) freelist.remove(result);
    }
    NodeType* findMin(){return head[0];}
    NodeType* findMax()
    {
        Node *result = 0, **tower = head;
        for(int level = currentLevel; level >= 0; --level)
            for(Node* node; node = tower[level]; tower = node->tower)
                result = node;
        return result;
    }
    class Iterator
    {
        Node* current;
    public:
        Iterator(Node* node): current(node){}
        Iterator& operator++()
        {
            assert(current);
            current = current->tower[0];
            return *this;
        }
        NodeType& operator*()
        {
            assert(current);
            return *current;
        }
        NodeType* operator->()
        {
            assert(current);
            return current;
        }
        bool operator==(Iterator const& rhs){return current == rhs.current;}
        bool operator!=(Iterator const& rhs){return current != rhs.current;}
    };
    Iterator begin(){return Iterator(findMin());}
    Iterator end(){return Iterator(0);}
};

}//end namespace
#endif
