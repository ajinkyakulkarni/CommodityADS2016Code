#ifndef EMBT_H
#define EMBT_H
#include "EMVector.h"
namespace igmdk{

template<typename KEY, typename ITEM> class EMBPlusTree
{
    enum{NULL_IO_POINTER = -1};
    typedef KVPair<KEY, long long> Key;
    typedef KVPair<KEY, ITEM> Record;
    enum{B = 2048, M = 2 * (2 + B/2/sizeof(Key)),
        L = 2 * (1 + B/2/sizeof(Record))};
    struct Node
    {//last item contains only a pointer, the first M - 1 items contain data
        int size;
        Key next[M];
        Node(): size(1) {next[0].value = NULL_IO_POINTER;}
        int findChild(KEY const& key)
        {
            int i = 0;
            while(i < size - 1 && key >= next[i].key) ++i;
            return i;
        }
    };
    struct Leaf
    {
        int size;
        long long next;
        Record records[L];
        Leaf(): size(0), next(NULL_IO_POINTER) {}
        int inclusiveSuccessorRecord(KEY const& key)
        {
            int i = 0;
            while(i < size && key > records[i].key) ++i;
            return i;
        }
    };
    long long leafIndex(long long index){return -(index + 2);}
    long long root;
    EMVector<Node> nodes;
    EMVector<Leaf> leaves;
    void splitInternal(long long index, int child)
    {
        Node parent = nodes[index];
        long long childIndex = parent.next[child].value;
        Node left = nodes[childIndex], right;
        //copy middle item key into parent
        for(int i = parent.size++; i > child; --i)
            parent.next[i] = parent.next[i - 1];
        parent.next[child].key = left.next[M/2 - 1].key;
        parent.next[child + 1].value = nodes.getSize();
        //move items starting from middle into right
        right.size = M/2 + 1;
        for(int i = 0; i < right.size; ++i)
            right.next[i] = left.next[i + M/2 - 1];
        left.size = M/2;
        nodes.append(right);
        nodes.set(left, childIndex);
        nodes.set(parent, index);
    }
    void splitLeaf(long long index, int child)
    {
        Node parent = nodes[index];
        long long childIndex = parent.next[child].value;
        Leaf left = leaves[leafIndex(childIndex)], right;
        //copy middle item key into parent
        for(int i = parent.size++; i > child; --i)
            parent.next[i] = parent.next[i - 1];
        parent.next[child].key = left.records[L/2].key;
        parent.next[child + 1].value = leafIndex(leaves.getSize());
        //move items starting from middle into right
        left.size = right.size = L/2;
        for(int i = 0; i < right.size; ++i)
            right.records[i] = left.records[i + L/2];
        right.next = left.next;
        left.next = leaves.getSize();
        leaves.append(right);
        leaves.set(left, leafIndex(childIndex));
        nodes.set(parent, index);
    }
public:
    EMBPlusTree(string const& keyFilename, string const& recordFilename,
        long long storedRoot = NULL_IO_POINTER, int extraItemsKey = 0,
        int extraItemsRecord = 0): root(storedRoot),
        nodes(keyFilename, sizeof(Node), extraItemsKey),
        leaves(recordFilename, sizeof(Leaf), extraItemsRecord) {}
    long long getRoot(){return root;}

    long long findLeaf(KEY const& key)
    {
        long long current = root;
        while(current >= 0)
        {
            Node node = nodes[current];
            current = node.next[node.findChild(key)].value;
        }
        return current;
    }
    ITEM find(KEY const& key, bool& status)
    {
        status = true;
        long long current = findLeaf(key);
        if(current != NULL_IO_POINTER)
        {
            Leaf leaf = leaves[leafIndex(current)];
            int i = leaf.inclusiveSuccessorRecord(key);
            if(i < leaf.size && key == leaf.records[i].key)
                return leaf.records[i].value;
        }
        status = false;
    }
    bool shouldSplit(long long node)
    {
        return node < NULL_IO_POINTER ?
            leaves[leafIndex(node)].size == L : nodes[node].size == M;
    }
    void insert(KEY const& key, ITEM const& value)
    {
        if(root == NULL_IO_POINTER)
        {
            root = leafIndex(leaves.getSize());
            leaves.append(Leaf());
        }
        else if(shouldSplit(root))
        {//check if need to split the root
            Node newRoot;
            newRoot.next[0].value = root;
            bool wasLeaf = root < NULL_IO_POINTER;
            root = nodes.getSize();
            nodes.append(newRoot);
            wasLeaf ? splitLeaf(root, 0) : splitInternal(root, 0);
        }
        long long index = root;
        while(index > NULL_IO_POINTER)
        {
            Node node = nodes[index];
            int childI = node.findChild(key),child = node.next[childI].value;
            if(shouldSplit(child))
            {//split children on they way down if needed
                child < NULL_IO_POINTER ? splitLeaf(index, childI) :
                    splitInternal(index, childI);
                if(key > nodes[index].next[childI].key)
                    child = nodes[index].next[childI + 1].value;
            }
            index = child;
        }
        //insert the item into the leaf
        Leaf leaf = leaves[leafIndex(index)];
        int i = leaf.inclusiveSuccessorRecord(key);
        if(i < leaf.size && key == leaf.records[i].key)
            leaf.records[i].value = value;
        else
        {
            for(int j = leaf.size++; j > i; --j)
                leaf.records[j] = leaf.records[j - 1];
            leaf.records[i] = Record(key, value);
        }
        leaves.set(leaf, leafIndex(index));
    }
    void remove(KEY const& key)
    {
        long long current = findLeaf(key);
        if(current != NULL_IO_POINTER)
        {
            Leaf leaf = leaves[leafIndex(current)];
            int i = leaf.inclusiveSuccessorRecord(key);
            if(i < leaf.size && key == leaf.records[i].key)
            {
                --leaf.size;
                for(int j = i; j < leaf.size; ++j)
                    leaf.records[j] = leaf.records[j + 1];
            }
            leaves.set(leaf, leafIndex(leafIndex(current)));
        }
    }
};

}
#endif
