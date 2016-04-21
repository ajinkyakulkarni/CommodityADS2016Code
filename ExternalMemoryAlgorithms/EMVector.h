#ifndef EMVECTOR_H
#define EMVECTOR_H
#include "File.h"
#include "../Utils/Vector.h"
#include "../Utils/Stack.h"
#include "../Utils/Sort.h"
#include "../Heaps/Heap.h"
#include "../Utils/Queue.h"
#include "../MiscAlgs/Misc.h"
#include <cmath>
using namespace std;

namespace igmdk{

class SingleBlockBuffer
{
    BlockFile& blockFile;
    Vector<char> block;
    long long loadedBlock;
    bool changed;
    void flush()
    {
        if(changed) blockFile.writeBlock(loadedBlock, block);
        changed = false;
    }
    void loadBlock(long long blockId)
    {
        if(blockId != loadedBlock)
        {
            flush();
            blockFile.readBlock(blockId, block);
            loadedBlock = blockId;
        }
    }
public:
    SingleBlockBuffer(BlockFile& theBlockFile): changed(false),
        loadedBlock(-1), blockFile(theBlockFile),
        block(blockFile.getBlockSize(), 0) {}
    ~SingleBlockBuffer(){flush();}
    void get(char* data, long long blockId, int start, int n)
    {
        loadBlock(blockId);
        for(int i = 0; i < n; ++i) data[i] = block[start + i];
    }
    void set(char* data, long long blockId, int start, int n)
    {
        loadBlock(blockId);
        for(int i = 0; i < n; ++i) block[start + i] = data[i];
        changed = true;
    }
};

//have bug somewhere can't use with vector!
class LRUBlockBuffer
{
    typedef pair<Vector<char>, bool> ITEM;//bool is changed flag
    BlockFile& blockFile;
    LRUCache<ITEM, long long> blocks;
    ITEM* loadBlock(long long blockId)
    {
        ITEM* item = blocks.read(blockId);
        if(!item)
        {
            long long* evictee = blocks.evicteeOnWrite(blockId);
            if(evictee)
            {
                ITEM* evicteeItem = blocks.read(*evictee);
                if(evicteeItem->second) flushWork(evicteeItem, *evictee);
            }
            Vector<char> block(blockFile.getBlockSize(), 0);
            blockFile.readBlock(blockId, block);
            blocks.write(blockId, make_pair(block, false));
            item = blocks.read(blockId);
        }
        return item;
    }
    void flushWork(ITEM* item, long long blockId)
    {
        if(item->second)
        {
            blockFile.writeBlock(blockId, item->first);
            item->second = false;
        }
    }
public:
    LRUBlockBuffer(BlockFile& theBlockFile, int cacheSize = 5):
        blocks(cacheSize), blockFile(theBlockFile) {}
    ~LRUBlockBuffer()
    {
        for(LRUCache<ITEM, long long>::Iterator i = blocks.begin();
            i != blocks.end(); ++i) flushWork(&i->value, i->key);
    }
    void get(char* data, long long blockId, int start, int n)
    {
        ITEM* item = loadBlock(blockId);
        for(int i = 0; i < n; ++i) data[i] = item->first[start + i];
    }
    void set(char* data, long long blockId, int start, int n)
    {
        ITEM* item = loadBlock(blockId);
        for(int i = 0; i < n; ++i) item->first[start + i] = data[i];
        item->second = true;
    }
};

template<typename POD, typename BUFFER = SingleBlockBuffer> class EMVector
{
    BlockFile blockFile;
    long long size;
    int itemsPerBlock(){return blockFile.getBlockSize()/sizeof(POD);}
    long long block(long long i){return i/itemsPerBlock();}
    long long index(long long i){return i % itemsPerBlock();}
    BUFFER buffer;
public:
    long long getSize(){return size;}
    EMVector(string const&filename, int blockSize = 2048, int extraItems = 0)
        : buffer(blockFile), blockFile(filename, blockSize),
        size(blockFile.getSize() * itemsPerBlock() - extraItems)
        {assert(blockSize % sizeof(POD) == 0);}
    long long extraItems()
        {return blockFile.getSize() * itemsPerBlock() - size;}
    void append(POD const& item)
    {
        ++size;
        if(extraItems() < 0) blockFile.appendEmptyBlock();
        set(item, size - 1);
    }
    void set(POD const& item, long long i)
    {
        assert(i >= 0 && i < size);
        char* data = (char*)&item;
        buffer.set(data, block(i), index(i) * sizeof(POD), sizeof(POD));
    }
    POD operator[](long long i)
    {
        assert(i >= 0 && i < size);
        POD result;
        char* data = (char*)&result;
        buffer.get(data, block(i), index(i) * sizeof(POD), sizeof(POD));
        return result;
    }
    void removeLast()
    {
        assert(size > 0);
        --size;
    }

    friend void IOSort(EMVector& vector)
    {
        {
            long long C = sqrt(vector.getSize() * vector.itemsPerBlock()),
                Q = vector.getSize()/C, lastQSize = vector.getSize() % C;
            EMVector temp("IOSortTempFile.igmdk");
            typedef KVPair<POD, long long> HeapItem;
            Heap<HeapItem, KVComparator<POD, long long> > merger;
            Vector<pair<Queue<POD>, long long> > buffers(Q + 1);
            for(long long i = 0, k = 0; i < Q + 1; ++i)
            {
                long long n = i == Q ? lastQSize : C;
                if(n > 0)
                {//sort each block and write the result to a temp vector
                    Vector<POD> buffer;
                    for(long long j = 0; j < n; ++j)
                        buffer.append(vector[k++]);
                    quickSort(buffer.getArray(), 0, buffer.getSize() - 1);
                    for(long long j = 0; j < n; ++j) temp.append(buffer[j]);
                    //record number of unmerged items in the block
                    buffers[i].second = n - 1;
                    //put smallest item of each block on the heap
                    merger.insert(HeapItem(temp[i * n], i));
                }
            }
            for(long long i = 0; i < vector.getSize(); ++i)
            {//merge
                long long q = merger.getMin().value;
                vector.set(merger.deleteMin().key, i);
                bool bufferIsEmpty = buffers[q].first.isEmpty();
                if(!bufferIsEmpty || buffers[q].second > 0)
                {
                    if(bufferIsEmpty)
                    {//refill
                        long long j = 0, next = q * C +
                            (q == Q ? lastQSize : C) - buffers[q].second;
                        while(j < vector.itemsPerBlock() &&
                            buffers[q].second-- > 0)
                            buffers[q].first.push(temp[next + j++]);
                    }
                    merger.insert(HeapItem(buffers[q].first.pop(), q));
                }
            }
        }
        File::remove("IOSortTempFile.igmdk");
    }
};

}
#endif
