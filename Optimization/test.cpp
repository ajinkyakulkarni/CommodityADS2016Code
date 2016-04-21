#include "MetaHeuristics.h"
#include "SearchAlgorithms.h"
#include "../ComputationalGeometry/Point.h"
#include "../Utils/Vector.h"
#include "../MiscAlgs/Misc.h"
#include "../RandomNumberGeneration/Random.h"
#include "../HashTable/LinearProbingHashTable.h"
#include "../Graphs/Graph.h"
using namespace igmdk;

struct TSP2ChangeMove
{
    int i, j;
    Vector<Point2>& points;
    Vector<int>& current;
    Vector<int> best;
    double currentEval, bestScore;
    TSP2ChangeMove(Vector<Point2>& thePoints, Vector<int>& order)
    :   points(thePoints), i(0), j(0), currentEval(fullScore()), current(order), best(order), bestScore(currentEval)
    {assert(points.getSize() > 1);}
    double eval(int from, int to)const
    {
        EuclideanDistance<Point2>::Distance d;
        return d(points[current[from]], points[current[to]]);
    }
    double proposeMove()
    {
        i = GlobalRNG.mod(current.getSize()-1);
        j = GlobalRNG.mod(current.getSize()-1);
        if(j < i) swap(i, j);
        double r = improvement();
        return r;
    }
    double currentScore(){return currentEval;}
    double bestEval(){return bestScore;}
    double nextScore()
    {
        double iFactor = i ? (eval(i-1,i)- eval(i-1,j)) : 0;
        double x = currentEval-iFactor-eval(j,j+1)+eval(i,j+1);
        return  x;
    }
    double improvement(){return currentEval - nextScore();}
    double fullScore()const
    {
        double result = 0;
        for(int i = 1; i < points.getSize(); ++i) result += eval(i-1, i);
        return result;
    }
    void applyMove()
    {
        currentEval = nextScore();
        current.reverse(i, j);
    }
    void updateBest()
    {
        if(currentEval < bestScore)
        {
            bestScore = currentEval;
            best = current;
        }
    }
    void bigMove()
    {
        updateBest();
        GlobalRNG.randomPermutation(current.getArray(), current.getSize());
        currentEval = fullScore();
    }
    void localSearchBest(){localSearch(*this, 1 << 10, 100);}
    bool isPruned(Vector<int> const& permutation)
    {
        //should use cost of points up i and cost of MST of remaining points
        return false;
    }
    void updateBest(Vector<int> const& permutation)
    {
        for(int i = 0; i < permutation.getSize(); ++i) DEBUG(permutation[i]);
        current = permutation;
        currentEval = fullScore();
        updateBest();
    }
};

struct TSPRandomInstance
{
    Vector<Point2> points;
    Vector<int> current;
    TSPRandomInstance(int N): current(N, 0)
    {
        for(int i = 0; i < N; ++i) current[i] = i;
        //to verify optimizality can just that order of points is the
        //same on the tour as on the point set's convex hull
        for(int i  = 0; i < N; ++i)
        {
            points.append(Point2(GlobalRNG.uniform01(), GlobalRNG.uniform01()));
        }
    }
};

void testSA2()
{
    int N = 10000;
    TSPRandomInstance instance(N);
    long long maxTries = 1<<20;
    TSP2ChangeMove move(instance.points, instance.current);
    double T = 10000, coolingFactor = 0.9999;
    simulatedAnnealing(move, T, coolingFactor, maxTries);
    DEBUG(move.currentScore());
}

void testLocalSearch2()
{
    int N = 10000;
    TSPRandomInstance instance(N);
    long long maxTries = 1ull<<20;
    TSP2ChangeMove move(instance.points, instance.current);
    localSearch(move, maxTries, 100);
    DEBUG(move.currentScore());
}

template<typename PROBLEM> class BranchAndBoundPermutation
{
    Vector<int> permutation;
    PROBLEM& problem;
    int N;
public:
    BranchAndBoundPermutation(int theN, PROBLEM& theProblem)
    :   N(theN), problem(theProblem){permutation.append(0);}
    typedef int Move;
    bool processSolution()
    {
        if(permutation.getSize() >= N)
        {
            problem.updateBest(permutation);
            return true;
        }
        return problem.isPruned(permutation);
    }
    Vector<KVPair<double, int> > generateMoves()
    {
        Vector<bool> isIncluded(N, false);
        for(int i = 0; i < permutation.getSize(); ++i)
            isIncluded[permutation[i]] = true;
        Vector<KVPair<double, int> > result;
        for(int i = 0; i < N; ++i)
            if(!isIncluded[i]) result.append(KVPair<double, int>(0, i));
        return result;
    }
    void move(int next){permutation.append(next);}
    void undoMove(int next){permutation.removeLast();}
};

void testBranchAndBound()
{
    int N = 4;
    TSPRandomInstance instance(N);
    TSP2ChangeMove problem(instance.points, instance.current);
    BranchAndBoundPermutation<TSP2ChangeMove> bb(N, problem);
    branchAndBound(bb);
    DEBUG(problem.bestScore);
}

void testIteratedLocalSearch2()
{
    int N = 10000;
    TSPRandomInstance instance(N);
    long long maxBigMoves = 1<<6;
    TSP2ChangeMove move(instance.points, instance.current);
    iteratedLocalSearch(move, maxBigMoves);
    DEBUG(move.bestEval());
}

struct GraphProblem
{
    GraphAA<double> graph;
    int from, to;
    int start(){return from;}
    bool isGoal(int i, LinearProbingHashTable<int, int>){return i == to;}
    Vector<int> nextStates(int j, LinearProbingHashTable<int, int>)
    {
        Vector<int> result;
        for(GraphAA<double>::AdjacencyIterator i = graph.begin(j); i != graph.end(j); ++i)
        {
            result.append(i.to());
        }
        return result;
    }
    double lowerbound(int i, LinearProbingHashTable<int, int>){return 0;}
    double distance(int k, int j)
    {
        for(GraphAA<double>::AdjacencyIterator i = graph.begin(j); i != graph.end(j); ++i)
        {
            if(i.to() == k) return i.data();
        }
        return 0;
    }
};

void timeSRT2()
{
    typedef GraphAA<double> G;
    G sp;
    GraphProblem Gp;
	for(int i = 0; i < 6; ++i)
	{
		Gp.graph.addVertex();
	}
	Gp.graph.addEdge(0,1,6);
	Gp.graph.addEdge(0,2,8);
	Gp.graph.addEdge(0,3,18);
	Gp.graph.addEdge(1,4,11);
	Gp.graph.addEdge(2,3,9);
	Gp.graph.addEdge(4,5,3);
	Gp.graph.addEdge(5,2,7);
	Gp.graph.addEdge(5,3,4);
	Gp.from = 0;
	Gp.to = 5;

	AStar<GraphProblem, int> dk(Gp);
	for(LinearProbingHashTable<int, int>::Iterator iter = dk.pred.begin();
        iter != dk.pred.end(); ++iter)
    {
        DEBUG(iter->key);
        DEBUG(iter->value);
    }
}

struct GraphProblem2
{
    GraphAA<double> graph;
    int from, to;
    int start(){return from;}
    bool isGoal(int i, Stack<int>){return i == to;}
    Vector<int> nextStates(int j, Stack<int>)
    {
        Vector<int> result;
        for(GraphAA<double>::AdjacencyIterator i = graph.begin(j); i != graph.end(j); ++i)
        {
            result.append(i.to());
        }
        return result;
    }
    double lowerbound(int i, Stack<int>){return 0;}
    double distance(int k, int j)
    {
        for(GraphAA<double>::AdjacencyIterator i = graph.begin(j); i != graph.end(j); ++i)
        {
            if(i.to() == k) return i.data();
        }
        return 0;
    }
};

void timeSRT3()
{
    typedef GraphAA<double> G;
    G sp;
    GraphProblem2 Gp;
	for(int i = 0; i < 6; ++i)
	{
		Gp.graph.addVertex();
	}
	Gp.graph.addEdge(0,1,6);
	Gp.graph.addEdge(0,2,8);
	Gp.graph.addEdge(0,3,18);
	Gp.graph.addEdge(1,4,11);
	Gp.graph.addEdge(2,3,9);
	Gp.graph.addEdge(4,5,3);
	Gp.graph.addEdge(5,2,7);
	Gp.graph.addEdge(5,3,4);
	Gp.from = 0;
	Gp.to = 5;

	RecursiveBestFirstSearch<GraphProblem2, int> dk(Gp);
	while(!dk.pred.isEmpty())
    {
        DEBUG(dk.pred.pop());
    }
}

int main()
{
    timeSRT2();
    timeSRT3();
    testSA2();
    testLocalSearch2();
    testIteratedLocalSearch2();
    testBranchAndBound();
	return 0;
}
