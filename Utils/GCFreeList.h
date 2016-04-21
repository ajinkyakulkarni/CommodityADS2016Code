#ifndef GC_FREE_LIST_H
#define GC_FREE_LIST_H
#include "Utils.h"
namespace igmdk{

template<typename ITEM> struct SimpleDoublyLinkedList
{
    struct Node
    {
        ITEM item;
        Node *next, *prev;
        template<typename ARGUMENT>
        Node(ARGUMENT const& argument, Node* theNext): item(argument),
            next(theNext), prev(0) {}
    } *root, *last;
    SimpleDoublyLinkedList(): root(0), last(0){}
    void prepend(Node* node)
    {
        assert(node);
        node->next = root;
        if(root) root->prev = node;
        node->prev = 0;
        root = node;
        if(!last) last = node;
    }
    void cut(Node* node)
    {
        assert(node);
        (node == last ? last : node->next->prev) = node->prev;
        (node == root ? root : node->prev->next) = node->next;
    }
    bool isEmpty(){return !root;}
    ~SimpleDoublyLinkedList()
    {
        while(root)
        {
            Node* toBeDeleted = root;
            root = root->next;
            delete toBeDeleted;
        }
    }
};

template<typename ITEM> struct StaticFreelist
{
    int capacity, size, maxSize;
    struct Item
    {
        ITEM item;
        union
        {
            Item* next;
            typename SimpleDoublyLinkedList<StaticFreelist>::Node* cameFrom;
        };
    } *nodes, *returned;
    StaticFreelist(int fixedSize): capacity(fixedSize), size(0), maxSize(0),
        returned(0), nodes(rawMemory<Item>(fixedSize)){}
    bool isFull(){return size >= capacity && !returned;}
    bool isEmpty(){return size <= 0;}
    Item* allocate()
    {
        Item* result = returned;
        if(result) returned = returned->next;
        else result = &nodes[maxSize++];
        ++size;
        return result;
    }
    void remove(Item* item)
    {
        item->item.~ITEM();
        item->next = returned;
        returned = item;
        --size;
    }
    ~StaticFreelist()
    {//O(1) if all items are returned, O(maxSize) otherwise
        if(!isEmpty())
        {//mark allocated nodes, unmark returned ones, destruct marked ones
            bool* toDelete = new bool[maxSize];
            for(int i = 0; i < maxSize; ++i) toDelete[i] = true;
            while(returned)
            {//nodes must come from this list for this to work
                toDelete[returned - nodes] = false;
                returned = returned->next;
            }
            for(int i = 0; i < maxSize; ++i)
                if(toDelete[i])nodes[i].item.~ITEM();
            delete[] toDelete;
        }
        rawDelete(nodes);
    }
};

template<typename ITEM> class Freelist
{
    enum{MAX_BLOCK_SIZE = 8192, MIN_BLOCK_SIZE = 8, DEFAULT_SIZE = 32};
    int totalSize;
    typedef SimpleDoublyLinkedList<StaticFreelist<ITEM> > ListType;
    typedef typename StaticFreelist<ITEM>::Item Item;
    typedef typename ListType::Node NodeType;
    ListType next, full;
    Freelist(Freelist const&);
    Freelist& operator=(Freelist const&);
public:
    Freelist(int initialSize = DEFAULT_SIZE): totalSize(max<int>(
        MIN_BLOCK_SIZE, min<int>(initialSize, MAX_BLOCK_SIZE))) {}
    ITEM* allocate()
    {
        if(next.isEmpty())
        {
            next.prepend(new NodeType(totalSize, 0));
            totalSize = min<int>(totalSize * 2, MAX_BLOCK_SIZE);
        }
        NodeType* root = next.root;
        Item* result = root->item.allocate();
        result->cameFrom = root;
        if(root->item.isFull())
        {
            next.cut(root);
            full.prepend(root);
        }
        return (ITEM*)result;
    }
    void remove(ITEM* item)
    {
        if(!item) return;
        Item* node = (Item*)(item);
        NodeType* cameFrom = node->cameFrom;
        StaticFreelist<ITEM>& bucket = cameFrom->item;
        bool wasFull = bucket.isFull();
        bucket.remove(node);
        if(bucket.isEmpty())//if 1 item buckets were allowed this would fail
        {//as item would be in full, not next
            totalSize -= bucket.capacity;
            next.cut(cameFrom);
            delete cameFrom;
        }
        else if(wasFull)
        {
            full.cut(cameFrom);
            next.prepend(cameFrom);
        }
    }
};
}//end namespace
#endif
