#ifndef MACHINELEARNINGOTHER_H
#define MACHINELEARNINGOTHER_H
#include "LearningCommon.h"
#include "../Utils/Utils.h"
#include "../MiscAlgs/Misc.h"
#include "../HashTable/LinearProbingHashTable.h"
#include "../ComputationalGeometry/KDTree.h"
#include "../ComputationalGeometry/Point.h"
#include "../NumericalMethods/Matrix.h"
#include "../RandomTreap/LCPTreap.h"
#include "../RandomNumberGeneration/Statistics.h"
#include <cmath>
namespace igmdk{

Matrix<double> randomProjection(int fromD, int toD)
{
    Matrix<double> result(toD, fromD);
    for(int i = 0; i < result.rows; ++i)
        for(int j = 0; j < result.columns; ++j)
            result(i, j) = GlobalRNG.bernoulli(0.5) ? -1 : 1;
    return result;
}

template<typename POINT, typename TREE = VpTree<POINT, int,
    typename EuclideanDistance<POINT>::Distance> > struct KMeans
{
    static Vector<int> findClusters(Vector<POINT>& points, int k,
        int maxIterations = 1000)
    {
        assert(k > 0 && k <= points.getSize() && points.getSize() > 0);
        //generate initial assignment
        Vector<int> assignments;
        //each cluster has at least 1 point, rest are random
        for(int i = 0; i < points.getSize(); ++i)
            assignments.append(i < k ? i : GlobalRNG.next() % k);
        bool converged = false;
        for(int m = 0; !converged && m < maxIterations; ++m)
        {//calculate centroids
            Vector<int> counts(k, 0);
            Vector<POINT> centroids(k, points[0] * 0);
            for(int i = 0; i < points.getSize(); ++i)
            {
                ++counts[assignments[i]];
                centroids[assignments[i]] += points[i];
            }
            for(int i = 0; i < k; ++i) centroids[i] *= 1.0/counts[i];
            TREE t;
            for(int i = 0; i < k; ++i) t.insert(centroids[i], i);
            //assign each point to the closest centroid
            converged = true;
            for(int i = 0; i < points.getSize(); ++i)
            {
                int best = t.nearestNeighbor(points[i])->value;
                if(best != assignments[i])
                {
                    converged = false;
                    assignments[i] = best;
                }
            }
        }
        return assignments;
    }
};

double UCB1(double averageValue, int nTries, int totalTries)
    {return averageValue + sqrt(2 * log(totalTries)/nTries);}

template<typename PROBLEM> void TDLearning(PROBLEM& p)
{
    while(p.hasMoreEpisodes())
    {
        double valueCurrent = p.startEpisode();
        while(!p.isInFinalState())
        {
            double valueNext = p.pickNextState();
            p.updateCurrentStateValue(p.learningRate() * (p.reward() +
                p.discountRate() * valueNext - valueCurrent));
            p.goToNextState();
            valueCurrent = valueNext;
        }
        p.updateCurrentStateValue(p.learningRate() *
            (p.reward() - valueCurrent));
    }
}

struct DiscreteValueFunction
{
    Vector<pair<double, int> > values;
    double learningRate(int state){return 1.0/values[state].second;}
    void updateValue(int state, double delta)
    {
        ++values[state].second;
        values[state].first += delta;
    }
    DiscreteValueFunction(int n): values(n, make_pair(0.0, 1)){}
};

struct LinearCombinationValueFunction
{
    Vector<double> weights;
    int n;
    double learningRate(){return 1.0/n;}
    void updateWeights(Vector<double> const& stateFeatures, double delta)
    {//set one of the state features to 1 to have a bias weight
        assert(stateFeatures.getSize() == weights.getSize());
        for(int i = 0; i < weights.getSize(); ++i)
            weights[i] += delta * stateFeatures[i];
        ++n;
    }
    LinearCombinationValueFunction(int theN): weights(theN, 0), n(1) {}
};

struct APriori
{
    LcpTreap<Vector<int>, int> counts;
    int processBasket(Vector<int> const& basket, int round,
        int rPrevMinCount = 0, int r1MinCount = 0)
    {
        int addedCount = 0;
        if(basket.getSize() > round)
        {
            Combinator c(round, basket.getSize());
            do//prepare the current combination of ids, needn't sort if each
            {//basket is already sorted
                Vector<int> key, single;
                for(int i = 0; i < round; ++i) key.append(basket[c.c[i]]);
                quickSort(key.getArray(), 0, key.getSize() - 1);
                int* count = counts.find(key);
                if(count) ++*count;//combination is frequent if already
                else if(round == 1)//frequent or round is 1
                {
                    counts.insert(key, 1);
                    ++addedCount;
                }
                else//combination is frequent if the last item and
                {//combination without the last item are both frequent
                    single.append(key.lastItem());
                    if(*counts.find(single) >= r1MinCount)
                    {
                        key.removeLast();
                        if(*counts.find(key) >= rPrevMinCount)
                        {
                            key.append(single[0]);
                            counts.insert(key, 1);
                            ++addedCount;
                        }
                    }
                }
            }while(!c.next());
        }
        return addedCount;
    }
    void noCutProcess(Vector<Vector<int> >const& baskets, int nRounds)
    {
        for(int k = 1; k <= nRounds; ++k)
            for(int i = 0; i < baskets.getSize(); ++i)
                processBasket(baskets[i], k);
    }
};

}//end namespace
#endif

