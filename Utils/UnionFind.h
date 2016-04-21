#ifndef UNION_FIND_H
#define UNION_FIND_H

#include "../Utils/Utils.h"
#include "../Utils/Vector.h"
#include "../RandomTreap/Treap.h"

namespace igmdk{

class UnionFind
{
    Vector<int> parent;//parent or negated size of the tree
public:
    UnionFind(int size): parent(size, -1){}
    int find(int n)
        {return parent[n] < 0 ? n : (parent[n] = find(parent[n]));}
    void join(int i, int j)
    {
        int parentI = find(i), parentJ = find(j);
        if(parentI != parentJ)
        {//parent[parentI] and parent[parentJ] are negative sizes
            if(parent[parentI] > parent[parentJ]) swap(parentI, parentJ);
            parent[parentI] += parent[parentJ];
            parent[parentJ] = parentI;
        }
    }
    bool areEquivalent(int i, int j){find(i) == find(j);}
    int subsetSize(int i){return -parent[find(i)];}
    void increaseSize(int newSize)
        {while(parent.getSize() < newSize) parent.append(-1);}
};

class IntervalSetUnion
{
    Treap<int, bool> treap;
public:
    int find(int i){return treap.isEmpty() ? -1 : treap.successor(i)->key;}
    void merge(int i){return treap.remove(i);}
    void split(int i){return treap.insert(i, 0);}
};

}
#endif
