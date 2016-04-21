#ifndef PRIORITY_QUEUE_H
#define PRIORITY_QUEUE_H
#include "../Utils/Vector.h"
#include "../HashTable/ChainingHashTable.h"
namespace igmdk{

template<typename ITEM>
struct ReportDefault{void operator()(ITEM& item, int i){}};

template<typename ITEM, typename COMPARATOR = DefaultComparator<ITEM>,
    typename REPORTER = ReportDefault<ITEM> > class Heap
{
    REPORTER r;
    int getParent(int i){return (i - 1)/2;}
    int getLeftChild(int i){return 2 * i + 1;}
    void moveUp(int i)
    {
        ITEM temp = items[i];
        for(int parent; i > 0 && c.isLess(temp, items[parent =
             getParent(i)]); i = parent) r(items[i] = items[parent], i);
        r(items[i] = temp, i);
    }
    void moveDown(int i)
    {
        ITEM temp = items[i];
        for(int child; (child = getLeftChild(i)) < items.getSize();
            i = child)
        {//find smaller child
            int rightChild = child + 1;
            if(rightChild < items.getSize() && c.isLess(items
                [rightChild], items[child])) child = rightChild;
            //replace with the smaller child if any
            if(!c.isLess(items[child], temp)) break;
            r(items[i] = items[child], i);
        }
        r(items[i] = temp, i);
    }
public:
    COMPARATOR c;
    Vector<ITEM> items;
    Heap(COMPARATOR const& theComparator = COMPARATOR(), REPORTER const&
        theReporter = REPORTER()): r(theReporter), c(theComparator) {}
    bool isEmpty(){return items.getSize() <= 0;}
    int getSize(){return items.getSize();}
    ITEM const& getMin()
    {
        assert(!isEmpty());
        return items[0];
    }
    void insert(ITEM const& item)
    {
        items.append(item);
        moveUp(items.getSize() - 1);
    }
    void changeKey(int i, ITEM const& item)
    {
        assert(i >= 0 && i <= items.getSize());
        bool decrease = c.isLess(item, items[i]);
        items[i] = item;
        decrease ? moveUp(i) : moveDown(i);
    }
    ITEM deleteMin(){return remove(0);}
    ITEM remove(int i)
    {
        assert(i >= 0 && i <= items.getSize());
        ITEM result = items[i];
        r(result, -1);
        if(items.getSize() > i)
        {
            items[i] = items.lastItem();
            r(items[i], i);
            moveDown(i);
        }
        items.removeLast();
        return result;
    }
};

template<typename ITEM, typename COMPARATOR = DefaultComparator<ITEM>,
    typename Handle = int> class IndexedHeap
{
    ChainingHashTable<Handle, int> map;
    typedef typename ChainingHashTable<int, int>::NodeType* POINTER;
    typedef KVPair<ITEM, POINTER> Item;
    typedef KVComparator<ITEM, POINTER, COMPARATOR> Comparator;
    struct Reporter
        {void operator()(Item& item, int i){item.value->value = i;}};
    Heap<Item, Comparator, Reporter> h;
public:
    IndexedHeap(COMPARATOR const& theComparator = COMPARATOR()):
        h(Comparator(theComparator)) {}
    int getSize(){return h.getSize();}
    ITEM* find(Handle handle)
    {
        int* pointer = map.find(handle);
        return pointer ? &h.items[*pointer].key : 0;
    }
    bool isEmpty(){return h.isEmpty();}
    void insert(ITEM const& item, Handle handle)
    {
        POINTER p = map.insert(handle, h.getSize());
        h.insert(Item(item, p));
    }
    ITEM const& getMin(){return h.getMin().key;}
    ITEM deleteMin()
    {
        Item result = h.deleteMin();
        map.remove(result.value->key);
        return result.key;
    }
    void changeKey(ITEM const& item, Handle handle)
    {
        POINTER p = map.findNode(handle);
        if(p) h.changeKey(p->value, Item(item, p));
        else insert(item, handle);
    }
    void deleteKey(Handle handle)
    {
        assert(find(handle));
        h.remove(*map.find(handle));
        map.remove(handle);
    }
};

template<typename ITEM, typename COMPARATOR = DefaultComparator<ITEM> >
class IndexedArrayHeap
{
    Vector<int> map;
    typedef KVPair<ITEM, int> Item;
    typedef KVComparator<ITEM, int, COMPARATOR> Comparator;
    struct Reporter
    {
        Vector<int>& pmap;
        Reporter(Vector<int>& theMap): pmap(theMap) {}
        void operator()(Item& item, int i){pmap[item.value] = i;}
    };
    Heap<Item, Comparator, Reporter> h;
public:
    IndexedArrayHeap(COMPARATOR const& theComparator = COMPARATOR()):
        h(Comparator(theComparator), Reporter(map)) {}
    ITEM* find(int handle)
    {
        int pointer = map[handle];
        return pointer != -1 ? &h.items[pointer].key : 0;
    }
    bool isEmpty(){return h.isEmpty();}
    void insert(ITEM const& item, int handle)
    {
        if(handle >= map.getSize())
            for(int i = map.getSize(); i <= handle; ++i) map.append(-1);
        h.insert(Item(item, handle));
    }
    ITEM const& getMin(){return h.getMin().key;}
    ITEM deleteMin()
    {
        Item result = h.deleteMin();
        map[result.value] = -1;
        return result.key;
    }
    void changeKey(ITEM const& item, int handle)
    {
        int p = map[handle];
        if(p != -1) h.changeKey(p, Item(item, handle));
        else insert(item, handle);
    }
    void deleteKey(int handle)
    {
        h.remove(map[handle]);
        map[handle] = -1;
    }
};
}
#endif
