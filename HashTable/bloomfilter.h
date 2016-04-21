#ifndef BLOOM_FILTER_H
#define BLOOM_FILTER_H

#include <cassert>

#include "../Utils/Bitset.h"
#include "HashFunction.h"
using namespace std;

namespace igmdk{

template<typename ITEM, typename HASHER = EHash32<MHash<PrimeHash> > >
class BloomFilter
{
    Bitset<unsigned char> items;
    HASHER h1, h2;
    int nHashes;
    int hash(int hash1, int hash2, int i)
    {
        if(i == 0) return hash1;
        if(i == 1) return hash2;
        return (hash1 + i * hash2) % items.getSize();
    }
public:
    BloomFilter(int size, int theNHashes = 7): nHashes(theNHashes), items(
        size), h1(size), h2(size) {assert(size > 0 && theNHashes > 0);}
    void insert(ITEM const& item)
    {
        int hash1 = h1.hash(item), hash2 = h2.hash(item);
        for(int i = 0; i < nHashes; ++i) items.set(hash(hash1, hash2, i));
    }
    bool isInserted(ITEM const& item)
    {
        int hash1 = h1.hash(item), hash2 = h2.hash(item);
        for(int i = 0; i < nHashes; ++i)
            if(!items[hash(hash1, hash2, i)]) return false;
        return true;
    }
};

}
#endif
