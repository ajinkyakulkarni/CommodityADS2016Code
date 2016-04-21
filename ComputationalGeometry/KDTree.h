#ifndef KDTREE_H
#define KDTREE_H
#include "../Utils/Utils.h"
#include "../Utils/Debug.h"
#include "../Utils/Vector.h"
#include "../Heaps/Heap.h"
#include "Point.h"
#include "../HashTable/ChainingHashTable.h"
#include "../RandomNumberGeneration/Random.h"
#include "../RandomNumberGeneration/Statistics.h"
#include "../Utils/GCFreelist.h"
#include <cmath>
namespace igmdk{

template<typename NODE> struct QNode
{
    NODE* node;
    double d;
    bool operator<(QNode const& rhs)const{return d > rhs.d;}
    static double dh(Heap<QNode>& heap, int k)
    {
        return heap.getSize() < k ?
            numeric_limits<double>::max() : heap.getMin().d;
    }
};

template<typename KEY, typename VALUE, int D,
    typename INDEXED_COMPARATOR = LexicographicComparator<KEY> > class KDTree
{
    INDEXED_COMPARATOR c;
    struct Node
    {
        KEY key;
        VALUE value;
        Node *left, *right;
        Node(KEY const& theKey, VALUE const& theValue): key(theKey),
            value(theValue), left(0), right(0) {}
    }* root;
    Freelist<Node> freelist;
public:
    typedef Node NodeType;
    bool isEmpty(){return !root;}
    KDTree(INDEXED_COMPARATOR theComparator = INDEXED_COMPARATOR()): root(0),
        c(theComparator) {}
    Node* constructFrom(Node* node)
    {
        Node* tree = 0;
        if(node)
        {
            tree = new(freelist.allocate())Node(node->key, node->value);
            tree->left = constructFrom(node->left);
            tree->right = constructFrom(node->right);
        }
        return tree;
    }
    KDTree(KDTree const& other): c(other.c)
        {root = constructFrom(other.root);}
    KDTree& operator=(KDTree const& rhs){return genericAssign(*this, rhs);}

    Node** findPointer(KEY const& key, Node*& parent)
    {
        Node* node, **pointer = &root;
        parent = 0;
        for(int i = 0; (node = *pointer) && !c.isEqual(key, node->key);
            i = (i + 1) % D)
        {
            parent = node;
            pointer = &(c.isLess(key, node->key, i) ?
                node->left : node->right);
        }
        return pointer;
    }
    VALUE* find(KEY const& key)
    {
        Node *node = *findPointer(key, node);
        return node ? &node->value : 0;
    }
    void insert(KEY const& key, VALUE const& value)
    {
        Node *dummy, **pointer = findPointer(key, dummy);
        if(*pointer) (*pointer)->value = value;
        else *pointer = new(freelist.allocate())Node(key, value);
    }

    void rangeQuery(KEY const& l, KEY const& u, bool* dimensions,
        Vector<Node*>& result, Node* node, int i)
    {
        if(!node) return;
        bool inRange = true;
        for(int j = 0; j < D; ++j)
            if(dimensions[j] && (c.isLess(node->key, l, j) ||
                c.isLess(u, node->key, i))) inRange = false;
        if(inRange) result.append(node);
        int j = (i + 1) % D;
        if(!(dimensions[i] && c.isLess(node->key, l, i)))
            rangeQuery(l, u, dimensions, result, node->left, j);
        if(!(dimensions[i] && c.isLess(u, node->key, i)))
            rangeQuery(l, u, dimensions, result, node->right, j);
    }
    void rangeQuery(KEY const& l, KEY const& u, bool* dimensions,
        Vector<Node*>& result)
        {rangeQuery(l, u, dimensions, result, root, 0);}

    template<typename DISTANCE> void distanceQuery(KEY const& x,
        double radius, Vector<Node*>& result, Node* node, int i,
        DISTANCE const& distance, bool cut)
    {
        if(!node) return;
        if(distance(node->key, x) <= radius) result.append(node);
        else if(cut) return;
        else cut = true;
        i = (i + 1) % D;
        Node* nodes[] = {node->left, node->right};
        for(int j = 0; j < 2; ++j)
            distanceQuery(x, radius, result, nodes[j], i, distance, cut);
    }
    template<typename DISTANCE> void distanceQuery(KEY const& x,
        double radius, Vector<Node*>& result, DISTANCE const& distance)
        {distanceQuery(x, radius, result, root, 0, distance, false);}

    typedef QNode<Node> HEAP_ITEM;
    template<typename DISTANCE> void kNN(Node* node, KEY const& key,
        Heap<HEAP_ITEM>& heap, int k, int i, KEY& partial,
        double partialDistance, DISTANCE const& distance)
    {
        double best = HEAP_ITEM::dh(heap, k);
        if(node && partialDistance < best)
        {//update partial distance
            double newPartialDistance = distance(key, node->key, i) -
                distance(key, partial, i);
            if(heap.getSize() < k)
            {
                HEAP_ITEM x = {node, distance(key, node->key)};
                heap.insert(x);
            }
            //use new partial distance to check for a cut again
            else if(newPartialDistance < best)
            {//incremental calculate-compare
                double d = distance(best, key, node->key);
                if(d < best)
                {
                    HEAP_ITEM x = {node, d};
                    heap.changeKey(0, x);
                }
            }
            int j = (i + 1) % D;
            //swap children for best order
            Node *l = node->left, *r = node->right;
            if(!c.isLess(key, node->key, i)) swap(l, r);
            kNN(l, key, heap, k, j, partial, partialDistance, distance);
            //set partial component to the node component, use the node
            //as temporary storage
            swap(partial[i], node->key[i]);
            kNN(r, key, heap, k, j, partial, newPartialDistance, distance);
            swap(partial[i], node->key[i]);
        }
    }
    template<typename DISTANCE> Vector<NodeType*> kNN(KEY const& key, int k,
        DISTANCE const& distance)
    {
        Heap<HEAP_ITEM> heap;
        KEY partial = key;
        kNN(root, key, heap, k, 0, partial, 0, distance);
        Vector<Node*> result;
        while(!heap.isEmpty()) result.append(heap.deleteMin().node);
        result.reverse();
        return result;
    }
    template<typename DISTANCE> NodeType* nearestNeighbor(KEY const& key,
        DISTANCE const& distance)
    {
        assert(!isEmpty());
        Node* parent, *result = *findPointer(key, parent);
        if(result) return result;
        Heap<HEAP_ITEM> heap;
        HEAP_ITEM x = {parent, distance(key, parent->key)};
        heap.insert(x);
        KEY partial = key;
        kNN(root, key, heap, 1, 0, partial, 0, distance);
        return heap.getMin().node;
    }
};

template<typename KEY, typename VALUE, typename DISTANCE> class VpTree
{
    DISTANCE distance;
    static double bound(double keyDistance, double rLow, double rHigh)
        {return max(0., max(keyDistance - rHigh, rLow - keyDistance));}
    struct Node
    {
        KEY key;
        VALUE value;
        double leftChildDistance, radius;
        Node *left, *right;
        Node(KEY const& theKey, VALUE const& theValue): key(theKey), left(0),
            right(0), value(theValue), leftChildDistance(0), radius(0) {}
        double leftChildBound(double keyDistance)
            {return bound(keyDistance, 0, leftChildDistance);}
        double rightChildBound(double keyDistance)
            {return bound(keyDistance, leftChildDistance, radius);}
    }* root;
    Freelist<Node> freelist;
public:
    typedef DISTANCE DISTANCE_TYPE;//update doc!
    DISTANCE const& getDistance(){return distance;}
    typedef Node NodeType;
    bool isEmpty(){return !root;}

    VpTree(DISTANCE const& theDistance = DISTANCE()): root(0),
        distance(theDistance) {}
    Node* constructFrom(Node* node)
    {
        Node* tree = 0;
        if(node)
        {
            tree = new(freelist.allocate())Node(node->key, node->value);
            tree->leftChildDistance = node->leftChildDistance;
            tree->radius = node->radius;
            tree->left = constructFrom(node->left);
            tree->right = constructFrom(node->right);
        }
        return tree;
    }
    VpTree(VpTree const& other): distance(other.distance)
        {root = constructFrom(other.root);}
    VpTree& operator=(VpTree const& rhs){return genericAssign(*this, rhs);}

    void distanceQuery(KEY const& key, double radius, Vector<Node*>& result,
        Node* node)
    {
        if(!node) return;
        double d = distance(node->key, key);
        if(d <= radius) result.append(node);
        if(node->leftChildBound(d) <= radius)
            distanceQuery(key, radius, result, node->left);
        if(node->rightChildBound(d) <= radius)
            distanceQuery(key, radius, result, node->right);
    }
    Vector<NodeType*> distanceQuery(KEY const& key, double radius)
    {
        Vector<NodeType*> result;
        distanceQuery(key, radius, result, root);
        return result;
    }

    typedef QNode<Node> HEAP_ITEM;
    void kNN(Node* node, KEY const& key, Heap<HEAP_ITEM>& heap, int k)
    {
        if(!node) return;
        //replace furthest node in heap with the current node if it's closer
        HEAP_ITEM x = {node, distance(key, node->key)};
        if(heap.getSize() < k) heap.insert(x);
        else if(x.d < HEAP_ITEM::dh(heap, k)) heap.changeKey(0, x);
        //expand closer child first
        double lb = node->leftChildBound(x.d),
            rb = node->rightChildBound(x.d);
        Node* l = node->left, *r = node->right;
        if(lb > rb)
        {
            swap(lb, rb);
            swap(l, r);
        }
        if(lb <= HEAP_ITEM::dh(heap, k)) kNN(l, key, heap, k);
        if(rb <= HEAP_ITEM::dh(heap, k)) kNN(r, key, heap, k);
    }
    Vector<NodeType*> kNN(KEY const& key, int k)
    {
        Heap<HEAP_ITEM> heap;
        kNN(root, key, heap, k);
        Vector<Node*> result;
        while(!heap.isEmpty()) result.append(heap.deleteMin().node);
        result.reverse();
        return result;
    }
    NodeType* nearestNeighbor(KEY const& key)
    {
        assert(!isEmpty());
        return kNN(key, 1)[0];
    }

    VALUE* find(KEY const& key)
    {
        Node* node = root;
        while(node && key != node->key)
            node = (distance(key, node->key) <= node->leftChildDistance) ?
                node->left : node->right;
        return node ? &node->value : 0;
    }
    void insert(KEY const& key, VALUE const& value)
    {
        Node **pointer = &root, *node;
        while((node = *pointer) && key != node->key)
        {
            double d = distance(key, node->key);
            node->radius = max(node->radius, d);
            if(!node->left) node->leftChildDistance = d;
            pointer = &(d <= node->leftChildDistance ?
                node->left : node->right);
        }
        if(node) node->value = value;
        else *pointer = new(freelist.allocate())Node(key, value);
    }
};

template<typename KEY, typename VALUE, typename DISTANCE> class KNNBruteForce
{
    DISTANCE distance;
    typedef KVPair<KEY, VALUE> Node;
    Vector<Node> nodes;
    struct QNode
    {
        double distance;
        int result;
        bool operator<(QNode const& rhs)const
            {return distance > rhs.distance;}
    };
public:
    KNNBruteForce(DISTANCE const& theDistance = DISTANCE()):
        distance(theDistance){}
    typedef Node NodeType;
    void insert(KEY const& key, VALUE const& value)
        {nodes.append(Node(key, value));}
    Vector<NodeType*> kNN(KEY const& key, int k)
    {
        Heap<QNode> q;
        for(int i = 0; i < nodes.getSize(); ++i)
        {
            QNode node = {distance(key, nodes[i].key), i};
            if(q.getSize() < k) q.insert(node);
            else if(node.distance < q.getMin().distance)
                q.changeKey(0, node);
        }
        Vector<NodeType*> result;
        while(!q.isEmpty()) result.append(&nodes[q.deleteMin().result]);
        result.reverse();
        return result;
    }
    NodeType* nearestNeighbor(KEY const& key){return kNN(key, 1)[0];}
};

class E2LSHHasher
{
    Vector<Xorshift64Hash> mappers;
    struct Hasher
    {
        Vector<double> a;
        double w, b;
        Hasher(int D, int r): w(r), b(GlobalRNG.uniform01() * w)
            {for(int i = 0; i < D; ++i) a.append(GlobalRNG.normal01());}
        int operator()(Vector<double> const& x)const
            {return int((a * x + b)/w);}
    };
    Vector<Hasher> h;
public:
    typedef unsigned long long RESULT_TYPE;
    typedef Vector<double> ITEM_TYPE;
    E2LSHHasher(int k, int l, int D, double w): mappers(l)
        {for(int i = 0; i < k * l; ++i) h.append(Hasher(D, w));}
    RESULT_TYPE operator()(ITEM_TYPE const& x, int bucket)const
    {
        Vector<int> result;
        int k = h.getSize()/mappers.getSize();
        for(int i = 0; i < k; ++i) {result.append(h[bucket * k + i](x)); //DEBUG(result.lastItem());
        }
        //system("PAUSE");
        return mappers[bucket].hash(result.getArray(), result.getSize());
    }
    static double p(double w, double r)
    {
        double z = r/w;
        return 2 * approxNormalCDF(z) - 1 -
            2/sqrt(2 * Random<>::PI())/z * (1 - exp(-z * z/2));
    }
    static double p1(double r){return p(r, r);}
    static double p2(double r, double c){return p(r, r * c);}
    static double distance(ITEM_TYPE const& x1, ITEM_TYPE const& x2)
    {
        EuclideanDistance<ITEM_TYPE>::Distance ed;
        return ed(x1, x2);
    }
};

namespace LSHKLFinder
{
    int LSHGetL(int k, double p1, double e)
    {
        double l = log(e)/log(1 - pow(p1, k));
        return (!isfinite(l) && l >= numeric_limits<int>::max()) ? -1 : 1 + int(l);
    }
    double LSHCost(int k, double e, double p1, double p2, int n)
    {
        int l = LSHGetL(k, p1, e);
        return (10 * k + pow(p2, k) * n) * l;
    }
    int minL(double p1, double e){return LSHGetL(1, p1, e);};
    int LSHFindK(double e, double p1, double p2, int n, int maxL)
    {
        int bestK = -1;
        double bestV;
        for(int k = 1;; ++k)
        {
            DEBUG(k);
            int l = LSHGetL(k, p1, e);
            DEBUG(l);
            double v = LSHCost(k, e, p1, p2, n);
            DEBUG(v);
            if(v < 0) break;

            if(bestK == -1 || (l > 0 && l < maxL && v < bestV)) {bestK = k; bestV = v;}
            if(l < 0 || l >= maxL) break;
        }
        DEBUG(bestV);
        //DEBUG(LSHGetL(bestK, p1, e));
        return bestK;
    }
}

template<typename HASHER> class LSH
{
    typedef typename HASHER::ITEM_TYPE ITEM;
    typedef typename HASHER::RESULT_TYPE RESULT_TYPE;
    Vector<ChainingHashTable<RESULT_TYPE, Vector<int> > > buckets;
    Vector<ITEM> items;
    HASHER g;
    double r2;
public:
    LSH(HASHER const& theG, int l, double theR2): buckets(l), g(theG), r2(theR2){}
    void insert(ITEM const& x)
    {
        for(int i = 0; i < buckets.getSize(); ++i)
        {
            typename HASHER::RESULT_TYPE hash = g(x, i);
            //DEBUG(i);
            //DEBUG(hash);
            Vector<int>* xBucket = buckets[i].find(hash);
            if(!xBucket)
            {
                buckets[i].insert(hash, Vector<int>());
                xBucket = buckets[i].find(hash);
            }
            xBucket->append(items.getSize());//have linear probing return chain instead?
        }
        items.append(x);
    }
    Vector<ITEM> cNeighbors(ITEM const& x)
    {
        Vector<ITEM> result;
        ChainingHashTable<int, bool> retrievedItems;
        int hitItems = 0;
        for(int i = 0; i < buckets.getSize(); ++i)
        {
            typename HASHER::RESULT_TYPE hash = g(x, i);
            //DEBUG(i);
            //DEBUG(hash);
            Vector<int>* xBucket = buckets[i].find(hash);
            if(xBucket)
                for(int i = 0; i < xBucket->getSize(); ++i)
                {
                    int itemIndex = (*xBucket)[i];
                    ++hitItems;
                    if(!retrievedItems.find(itemIndex))
                    {
                        retrievedItems.insert(itemIndex, true);
                        if(HASHER::distance(x, items[itemIndex]) < r2)
                            result.append(items[itemIndex]);
                    }
                }
        }
        //DEBUG(hitItems);
        return result;
    }
};

LSH<E2LSHHasher> buildE2LSH(int D, double r, double c, int maxL, double e = 10e-6, int maxN = 1000000)
{
    double p1 = E2LSHHasher::p(1, 1), r2 = r * (1 + c);
    int k = LSHKLFinder::LSHFindK(e, p1, E2LSHHasher::p(r, r2), maxN, maxL);
    //DEBUG(k);
    int l = LSHKLFinder::LSHGetL(k, p1, e);
    //DEBUG(l);
    return LSH<E2LSHHasher>(E2LSHHasher(k, l, D, r), l, r2);
}

template<typename HASHER> class NearestNeighborLSH
{
    typedef typename HASHER::ITEM_TYPE ITEM;
    Vector<LSH<HASHER> > lshs;//items are duplicated dont store them!
public:
    void addLSH(LSH<HASHER> const& lsh){lshs.append(lsh);}
    void insert(ITEM const& x)
        {for(int i = 0; i < lshs.getSize(); ++i) lshs[i].insert(x);}
    pair<ITEM, bool> cNeighbor(ITEM const& x)
    {
        for(int i = 0; i < lshs.getSize(); ++i)
        {
            Vector<ITEM> items = lshs[i].cNeighbors(x);
            if(items.getSize() > 0)
            {
                int best = -1, bestD;
                for(int j = 0; j < items.getSize(); ++j)
                {
                    double d = HASHER::distance(x, items[j]);
                    if(best == -1 || d < bestD)
                    {
                        best = j;
                        bestD = d;
                    }
                }
                return pair<ITEM, bool>(items[best], true);
            }
        }
        return pair<ITEM, bool>(ITEM(), false);
    }
};

NearestNeighborLSH<E2LSHHasher> buildE2NNLSH(int D, double rMin, double rMax, int maxL, double c = 1, double e = 10e-6, int maxN = 1000000)
{
    NearestNeighborLSH<E2LSHHasher> result;
    for(double r = rMin; r < rMax; r *= (1 + c))
    {
        result.addLSH(buildE2LSH(D, r, c, maxL, e, maxN));
    }
    return result;
}

}//end namespace
#endif
