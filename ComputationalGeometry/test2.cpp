#include <iostream>
#include <cmath>
#include "KDTree.h"
#include "../ComputationalGeometry/Point.h"
#include "../RandomNumberGeneration/Random.h"
#include "../NumericalMethods/NumericalMethods.h"
using namespace igmdk;

void testKDTree()
{
    KDTree<Point<int>, int, 2> kdtree;
    int N = 1500000;
    for(int i = 0; i < N; ++i)
    {
        int p1 = (GlobalRNG.next() % 1000), p2 = (GlobalRNG.next() % 1000);
        kdtree.insert(Point<int>(min(p1, p2), max(p1, p2)), i);
    }
    bool dimensions[2];
    dimensions[0] = true;
    dimensions[1] = true;
    for(int k = 0; k < 1; ++k)
    {
        Vector<KDTree<Point<int>, int, 2>::NodeType*> result;
        int point = 0;
        kdtree.rangeQuery(Point<int>(-999999999, point), Point<int>(point, 999999999), dimensions, result);
        DEBUG(result.getSize());
    }
}

void testKDTree2()
{
    KDTree<Point<int>, int, 2> kdtree;
    int N = 1500000;
    for(int i = 0; i < N; ++i)
    {
        kdtree.insert(Point<int>(GlobalRNG.next(), GlobalRNG.next()), i);
    }
    int M = 1500;
    for(int i = 0; i < M; ++i)
    {
        assert(kdtree.nearestNeighbor(Point<int>(GlobalRNG.next(), GlobalRNG.next()), EuclideanDistance<Point<int> >::DistanceIncremental()));
        int k = 2;
        assert(kdtree.kNN(Point<int>(GlobalRNG.next(), GlobalRNG.next()), k, EuclideanDistance<Point<int> >::DistanceIncremental()).getSize() == k);
    }
}

void testVpTree()
{
    VpTree<Point<int>,int, EuclideanDistance<Point<int> >::DistanceIncremental> tree;
    int N = 5;
    for(int i = 0; i < N; ++i)
    {
        tree.insert(Point<int>(i, i), i);
    }
    for(int i = 0; i < N; ++i)
    {
        int* result = tree.find(Point<int>(i, i));
        if(result) DEBUG(*result);
        assert(result && *result == i);
    }

    Vector<VpTree<Point<int>,int, EuclideanDistance<Point<int> >::DistanceIncremental>::NodeType*> result2 = tree.distanceQuery(Point<int>(0, 4), sqrt(10));

    for(int i = 0; i < result2.getSize(); ++i)
    {
        DEBUG(result2[i]->key[0]);
        DEBUG(result2[i]->key[1]);
    }

}

void testVpTree2()
{
    VpTree<Point<int>, int, EuclideanDistance<Point<int> >::DistanceIncremental> tree;
    int N = 1500000;
    for(int i = 0; i < N; ++i)
    {
        tree.insert(Point<int>(GlobalRNG.next(), GlobalRNG.next()), i);
    }
    int M = 15000;
    for(int i = 0; i < M; ++i)
    {
        assert(tree.nearestNeighbor(Point<int>(GlobalRNG.next(), GlobalRNG.next())));
        int k = 100;
        assert(tree.kNN(Point<int>(GlobalRNG.next(), GlobalRNG.next()), k).getSize() == k);
    }
}

struct FunctionTester
{
    void operator()()const
    {
        testVpTree2();
        testVpTree();
        testKDTree();
        testKDTree2();
    }
};

void testLSH()
{
    //DEBUG(1/exp(1));
    //DEBUG(E2LSHHasher::p(1, 1));
    int D = 2;
    LSH<E2LSHHasher> tree = buildE2LSH(D, 1, 1, D * 5);
    int N = 1000000;
    for(int i = 0; i < N; ++i)
    {
        tree.insert(Vector<double>(2, i));
    }
    for(int i = 0; i < 1; ++i)
    {
        Vector<Vector<double> > neighbors = tree.cNeighbors(Vector<double>(2, i));
        DEBUG(neighbors.getSize());
        for(int j = 0; j < neighbors.getSize(); ++j)
        {
            DEBUG(j);
            for(int k = 0; k < neighbors[j].getSize(); ++k) DEBUG(neighbors[j][k]);
        }
    }
}

void testLSH2()
{
    //DEBUG(1/exp(1));
    //DEBUG(E2LSHHasher::p(1, 1));
    int D = 2;
    NearestNeighborLSH<E2LSHHasher> tree = buildE2NNLSH(D, 1, 10, D * 5);
    int N = 100000;
    for(int i = 0; i < N; ++i)
    {
        tree.insert(Vector<double>(2, i));
    }
    int noneCount = 0;
    for(int i = 0; i < N; ++i)
    {
        pair<Vector<double>, bool> neighbor = tree.cNeighbor(Vector<double>(2, i));
        //DEBUG(neighbor.second);
        if(neighbor.second)
        {
            //for(int k = 0; k < neighbor.first.getSize(); ++k) DEBUG(neighbor.first[k]);
        }
        else ++noneCount;
    }
    DEBUG(noneCount);
}

void testLSH3()
{
    //DEBUG(1/exp(1));
    //DEBUG(E2LSHHasher::p(1, 1));
    int D = 100;
    NearestNeighborLSH<E2LSHHasher> tree = buildE2NNLSH(D, 1, 10, D * 0.5);
    int N = 100000;
    for(int i = 0; i < N; ++i)
    {
        Vector<double> x;
        for(int j = 0; j < D; ++j) x.append(i);
        tree.insert(x);
    }
    int noneCount = 0;
    for(int i = 0; i < N; ++i)
    {
        Vector<double> x;
        for(int j = 0; j < D; ++j) x.append(i);
        pair<Vector<double>, bool> neighbor = tree.cNeighbor(x);
        //DEBUG(neighbor.second);
        if(neighbor.second)
        {
            //for(int k = 0; k < neighbor.first.getSize(); ++k) DEBUG(neighbor.first[k]);
        }
        else ++noneCount;
    }//20 secs
    DEBUG(noneCount);
}

void testKD3()
{
    KDTree<Point<double, 100>, bool, 100> kdtree;
    int D = 100;
    int N = 100000;
    for(int i = 0; i < N; ++i)
    {
        Point<double, 100> x;
        for(int j = 0; j < D; ++j) x[j] = j;
        kdtree.insert(x, true);
    }
    for(int i = 0; i < N; ++i)
    {
        Point<double, 100> x;
        for(int j = 0; j < D; ++j) x[j] = j + GlobalRNG.uniform01();
        assert(kdtree.nearestNeighbor(x, EuclideanDistance<Point<double, 100> >::DistanceIncremental()));
    }
}
//LSH not better then KD-tree even for ML because index memory needed > data memory? Though much better than brute force!
//Before release or use must find data set or use case where it's much better!
void testKNNBF()
{
    KNNBruteForce<Point<double, 100>, bool, EuclideanDistance<Point<double, 100> >::Distance> kdtree;
    int D = 100;
    int N = 100000;
    for(int i = 0; i < N; ++i)
    {
        Point<double, 100> x;
        for(int j = 0; j < D; ++j) x[j] = j;
        kdtree.insert(x, true);
    }
    for(int i = 0; i < N; ++i)
    {
        Point<double, 100> x;
        for(int j = 0; j < D; ++j) x[j] = j + GlobalRNG.uniform01();
        assert(kdtree.nearestNeighbor(x));
    }
}

void DDDVPTree()
{
    Random<> rng(0);
    VpTree<Point<double>, int, EuclideanDistance<Point<double> >::DistanceIncremental> VPTree0to9;
    int D = 2;
    for(int i = 0; i < 10; ++i)
    {
        Point<double, 2> x;
        for(int j = 0; j < D; ++j) x[j] = rng.uniform01();
        VPTree0to9.insert(x, i);
    }

    cout << "breakpoint" << endl;
}

void DDDKDTree()
{
    Random<> rng(0);
    KDTree<Point<double, 2>, int, 2> KDTree0to9;
    int D = 2;
    for(int i = 0; i < 10; ++i)
    {
        Point<double, 2> x;
        for(int j = 0; j < D; ++j) x[j] = rng.uniform01();
        KDTree0to9.insert(x, i);
    }

    cout << "breakpoint" << endl;
}

int main()
{
    DDDVPTree();
    DDDKDTree();
    for(double i = 1; i <= 8; i++)
    {
        double p1 = E2LSHHasher::p(1, i);
        double p2 = E2LSHHasher::p(1.5, i);
        DEBUG(i);
        DEBUG(p1);
        DEBUG(p2);
        int l = 50;
        int k = 20;
        //for(int k = 1; k <= l; ++k)
        {
            //LSHKLFinder::LSHGetL(k, p1, 0.1);
            //if(l > 100) break;
            double p1Real = 1 - pow((1 - pow(p1, k)), l);
            double p2Real = 1 - pow((1 - pow(p2, k)), l);
            //if(p1Real < 0.9 || p2Real > 0.5) continue;
            DEBUG(k);
            DEBUG(l);
            DEBUG(p1Real);
            DEBUG(p2Real);
        }
        system("PAUSE");
    }
    //testLSH();
    //testLSH2();
    //testLSH3();
    //testKD3();//very fast
    //testKNNBF();//very slow
    return 0;
    /*NormalSummary result = MonteCarlo::simulate(SpeedTester<FunctionTester>(), 1);
    DEBUG(result.minimum);
    DEBUG(result.maximum);
    DEBUG(result.mean);
    DEBUG(result.variance);
    DEBUG(result.confidence9973());*/
	return 0;
}
