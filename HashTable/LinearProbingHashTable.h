#ifndef LINEAR_PROBING_HASH_TABLE_H
#define LINEAR_PROBING_HASH_TABLE_H
#include <new>
#include "HashFunction.h"
namespace igmdk{

template<typename KEY, typename VALUE, typename HASHER = EHash32<BUHash>,
typename COMPARATOR = DefaultComparator<KEY> > class LinearProbingHashTable
{
    int capacity, size;
    typedef KVPair<KEY, VALUE> Node;
    Node* table;
    bool* isOccupied;
    HASHER h;
    COMPARATOR c;
    void allocateTable(int requestedSize)
    {
        int bits = lgCeiling(max(requestedSize, 8));
        capacity = twoPower(bits);
        h = HASHER(bits);
        size = 0;
        table = rawMemory<Node>(capacity);
        isOccupied = new bool[capacity];
        for(int i = 0; i < capacity; ++i) isOccupied[i] = false;
    }
    static void cleanUp(Node* theTable, int theCapacity, bool* isOccupied)
    {
        for(int i = 0; i < theCapacity; ++i)
            if(isOccupied[i]) theTable[i].~Node();
        rawDelete(theTable);
        delete[] isOccupied;
    }
    void destroy(int cell)
    {
        table[cell].~Node();
        isOccupied[cell] = false;
        --size;
    }
    void resize()
    {
        int oldCapacity = capacity;
        Node* oldTable = table;
        bool* oldIsOccupied = isOccupied;
        allocateTable(2 * size);
        for(int i = 0; i < oldCapacity; ++i)
            if(oldIsOccupied[i]) insert(oldTable[i].key, oldTable[i].value);
        cleanUp(oldTable, oldCapacity, oldIsOccupied);
    }
public:
    typedef Node NodeType;
    void getSize(){return size;}
    LinearProbingHashTable(int initialCapacity = 8,
        COMPARATOR const& theComparator = COMPARATOR()): c(theComparator)
        {allocateTable(initialCapacity);}
    LinearProbingHashTable(LinearProbingHashTable const& rhs):
        capacity(rhs.capacity), h(rhs.h), size(rhs.size), c(rhs.c),
        isOccupied(new bool[capacity]), table(rawMemory<Node>(capacity))
    {
        for(int i = 0; i < capacity; ++i)
            if(isOccupied[i] = rhs.isOccupied[i]) table[i] = rhs.table[i];
    }
    LinearProbingHashTable& operator=(LinearProbingHashTable const& rhs)
        {return genericAssign(*this, rhs);}
    ~LinearProbingHashTable(){cleanUp(table, capacity, isOccupied);}
    int findNode(KEY const& key)
    {
        int cell = h.hash(key);
        for(;isOccupied[cell] && !c.isEqual(key, table[cell].key);
            cell = (cell + 1) % capacity);
        return cell;
    }
    VALUE* find(KEY const& key)
    {
        int cell = findNode(key);
        return isOccupied[cell] ? &table[cell].value : 0;
    }
    void insert(KEY const& key, VALUE const& value)
    {
        int cell = findNode(key);
        if(isOccupied[cell]) table[cell].value = value;
        else
        {
            new(&table[cell])Node(key, value);
            isOccupied[cell] = true;
            if(++size > capacity * 0.8) resize();
        }
    }
    void remove(KEY const& key)
    {
        int cell = findNode(key);
        if(isOccupied[cell])
        {//reinsert subsequent nodes in the found value's chain
            destroy(cell);
            while(isOccupied[cell = (cell + 1) % capacity])
            {
                Node temp = table[cell];
                destroy(cell);
                insert(temp.key, temp.value);
            }
            if(size < capacity * 0.1) resize();
        }
    }
    struct Iterator
    {
        int i;
        LinearProbingHashTable& t;
        void advance(){while(i < t.capacity && !t.isOccupied[i]) ++i;}
    public:
        Iterator(LinearProbingHashTable& theHashTable): i(0),
            t(theHashTable) {advance();}
        Iterator& operator++()
        {
            ++i;
            advance();
            return *this;
        }
        NodeType& operator*()const
            {assert(i < t.capacity); return t.table[i];}
        NodeType* operator->()const
            {assert(i < t.capacity); return &t.table[i];}
        bool operator!=(Iterator const& rhs)const{return i != rhs.i;}
    };
    Iterator begin(){return Iterator(*this);}
    Iterator end()
    {
        Iterator result(*this);
        result.i = capacity;
        return result;
    }
};

}
#endif
