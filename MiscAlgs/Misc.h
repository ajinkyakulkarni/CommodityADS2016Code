#ifndef MISC_H
#define MISC_H
#include "../Utils/Bitset.h"
#include "../Utils/Sort.h"
#include "../Utils/GCFreelist.h"
#include "../HashTable/LinearProbingHashTable.h"
#include <cmath>

namespace igmdk{

class CRC32
{
    unsigned int polynomial, constant[256];
public:
    CRC32(unsigned int thePolynomial = 0xFA567D89u):polynomial(thePolynomial)
    {
        for(int i = 0; i < 256; ++i)
        {
            constant[i] = i << 24;
            for(int j = 0; j < 8; ++j) constant[i] =
                (constant[i] << 1) ^ (constant[i] >> 31 ? polynomial : 0);
        }
    }
    unsigned int hash(unsigned char* array, int size, unsigned int crc = 0)
    {
        for(int i = 0; i < size; ++i)
            crc = (crc << 8) ^ constant[(crc >> 24) ^ array[i]];
        return crc;
    }
};

template<typename VALUE, typename KEY = int,
    typename HASHER = EHash32<BUHash> > class LRUCache
{
    typedef KVPair<KEY, VALUE> ITEM;
    typedef SimpleDoublyLinkedList<ITEM> LIST;
    typedef typename LIST::Node NODE;
    LIST l;
    int size, capacity;
    LinearProbingHashTable<KEY, NODE*, HASHER> h;
public:
    LRUCache(int theCapacity): size(0), capacity(theCapacity)
        {assert(capacity > 0);}
    struct Iterator
    {
        NODE* current;
        Iterator(NODE* node): current(node){}
        Iterator& operator++()
        {
            assert(current);
            current = current->next;
            return *this;
        };
        ITEM& operator*()const{assert(current); return current->item;}
        ITEM* operator->()const{assert(current); return &current->item;}
        bool operator!=(Iterator const& rhs)const
            {return current != rhs.current;}
    };
    VALUE* read(KEY const& k)
    {
        NODE** n = h.find(k);
        if(n)
        {
            assert(*n);
            l.cut(*n);
            l.prepend(*n);
            return &(*n)->item.value;
        }
        return 0;
    }
    KEY* evicteeOnWrite(KEY const& k)
    {
        if(size == capacity) return &l.last->item.key;
        return 0;
    }
    void write(KEY const& k, VALUE const& v)
    {
        KEY* evictee = evicteeOnWrite(k);
        NODE* n = l.last;
        if(evictee)
        {
            h.remove(*evictee);
            l.cut(n);
            n->item.key = k;
            n->item.value = v;
        }
        else
        {
            ++size;
            n = new NODE(ITEM(k, v), 0);
        }
        l.prepend(n);
        h.insert(k, l.root);
    }
    Iterator begin(){return Iterator(l.root);}
    Iterator end(){return Iterator(l.last);}
};

template<typename KEY, typename VALUE, typename HASHER, typename RESOURCE>
class ReadLRUCache
{
    RESOURCE const& r;
    LRUCache<KEY, VALUE, HASHER> c;
public:
    ReadLRUCache(RESOURCE const& theR, HASHER const& h, int capacity):
        r(theR), c(h, capacity) {}
    static VALUE* readWork(KEY const& k, LRUCache<KEY, VALUE, HASHER>& c,
        RESOURCE const& r)
    {
        VALUE* v = c.read(k);
        if(!v)
        {
            v = r.read(k);
            if(v) c.write(k, *v);
        }
        return v;
    }
    VALUE* read(KEY const& k){return readWork(k, c, r);}
};

template<typename KEY, typename VALUE, typename HASHER, typename RESOURCE>
class InstantCommitLRUCache
{
    RESOURCE & r;
    LRUCache<KEY, VALUE, HASHER> c;
public:
    InstantCommitLRUCache(RESOURCE& theR, HASHER const& h, int capacity):
        r(theR), c(h, capacity) {}
    VALUE* read(KEY const& k)
    {
        return ReadLRUCache<KEY, VALUE, HASHER, RESOURCE>::readWork(k, c, r);
    }
    void write(KEY const& k, VALUE const& v)
    {
        c.write(k, v);
        r.write(k, v);
    }
};

template<typename KEY, typename VALUE, typename HASHER, typename RESOURCE>
class DelayedCommitLRUCache
{
    typedef pair<VALUE, bool> MARKED_VALUE;
    RESOURCE & r;
    LRUCache<KEY, MARKED_VALUE, HASHER> c;
public:
    DelayedCommitLRUCache(RESOURCE& theR, HASHER const& h, int capacity):
        r(theR), c(h, capacity) {}
    VALUE* read(KEY const& k)
    {
        MARKED_VALUE* mv = ReadLRUCache<KEY, MARKED_VALUE, HASHER,
            RESOURCE>::readWork(k, c, r);
        return mv ? &mv->first : 0;
    }
    void write(KEY const& k, VALUE const& v)
    {
        MARKED_VALUE* mv = c.evicteeOnWrite(k);
        if(mv && mv->second) r.write(k, mv->first);
        c.write(k, MARKED_VALUE(v, true));
    }
    ~DelayedCommitLRUCache()
    {
        typedef typename LRUCache<KEY, MARKED_VALUE, HASHER>::Iterator I;
        for(I i = c.begin(); i != c.end(); ++i)
            if(i->second) r.write(i->first.key, i->first.value);
    }
};

class PrimeTable
{
    long long maxN;
    Bitset<> table;//marks odd numbers starting from 3
    long long nToI(long long n){return (n - 3)/2;}
public:
    PrimeTable(long long primesUpto): maxN(primesUpto - 1),
        table(nToI(maxN) + 1)
    {
        assert(primesUpto > 1);
        table.setAll(true);
        for(long long i = 3; i <= sqrt(maxN); i += 2)
            if(isPrime(i))//set every odd multiple i <= k <= maxN/i to false
                for(long long k = i; i * k <= maxN; k += 2)
                    table.set(nToI(i * k), false);
    }
    bool isPrime(long long n)
    {
        assert(n <= maxN);
        return n == 2 || (n > 2 && n % 2 && table[nToI(n)]);
    }
};

struct Permutator
{
    Vector<int> p;
    Permutator(int size){for(int i = 0; i < size; ++i) p.append(i);}
    bool next()
    {//find largest i such that p[i] < p[i + 1]
        int j = p.getSize() - 1, i = j - 1;
        while(i >= 0 && p[i] >= p[i + 1]) --i;
        bool backToIdentity = i == -1;
        if(!backToIdentity)
        {//find j such that p[j] is next largest element after p[i]
            while(i < j && p[i] >= p[j]) --j;
            swap(p[i], p[j]);
        }
        p.reverse(i + 1, p.getSize() - 1);
        return backToIdentity;//true if returned to smallest permutation
    }
    bool advance(int i)
    {
        assert(i >= 0 && i < p.getSize());
        quickSort(p.getArray(), i + 1, p.getSize() - 1,
            ReverseComparator<int>());
        return next();
    }
};

struct Combinator
{
    int n;
    Vector<int> c;
    Combinator(int m, int theN): n(theN), c(m, -1)
    {
        assert(m <= n && m > 0);
        skipAfter(0);
    }
    void skipAfter(int i)
    {//increment c[i] and reset all c[j] for j > i
        assert(i >= 0 && i < c.getSize());
        ++c[i];
        for(int j = i + 1; j < c.getSize(); ++j) c[j] = c[j - 1] + 1;
    }
    bool next()
    {//find rightmost c[i] which can be increased
        int i = c.getSize() - 1;
        while(i >= 0 && c[i] == n - c.getSize() + i) --i;
        bool finished = i == -1;
        if(!finished) skipAfter(i);
        return finished;
    }
};

struct Partitioner
{
    Vector<int> p;
    Partitioner(int n): p(n, 0) {assert(n > 0);}
    bool skipAfter(int k)
    {//set trailing elements to maximum values and call next
        assert(k >= 0 && k < p.getSize());
        for(int i = k; i < p.getSize(); ++i) p[i] = i;
        return next();
    }
    bool next()
    {//find rightmost p[j] which can be increased
        int m = 0, j = -1;
        for(int i = 0; i < p.getSize(); ++i)
        {
            if(p[i] < m) j = i;
            m = max(m, p[i] + 1);
        }
        bool finished = j == -1;
        if(!finished)
        {//increase it and reset the tail
            ++p[j];
            for(int i = j + 1; i < p.getSize(); ++i) p[i] = 0;
        }
        return finished;
    }
};

}
#endif
