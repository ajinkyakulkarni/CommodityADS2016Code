#ifndef SEARCH_ALGORITHMS_H
#define SEARCH_ALGORITHMS_H
#include "../Utils/Vector.h"
#include "../HashTable/LinearProbingHashTable.h"
#include "../Heaps/Heap.h"
#include "../Utils/Sort.h"
namespace igmdk{

template<typename PROBLEM, typename STATE_ID = unsigned long long>
struct AStar
{
    LinearProbingHashTable<STATE_ID, STATE_ID> pred;
    AStar(PROBLEM& p)
    {
        typedef KVPair<double, STATE_ID> QNode;
        IndexedHeap<QNode, KVComparator<double, STATE_ID> > pQ;
        STATE_ID j = p.start();
        pred.insert(j, -1);
        pQ.insert(QNode(p.lowerbound(j, pred), j), j);
        while(!pQ.isEmpty() && !p.isGoal(j = pQ.getMin().value, pred))
        {
            double dj = pQ.deleteMin().key - p.lowerbound(j, pred);
            Vector<STATE_ID> next = p.nextStates(j, pred);
            for(int i = 0; i < next.getSize(); ++i)
            {
                STATE_ID to = next[i];
                if(!pred.find(to)) pred.insert(to, j);
                double newDistance = dj + p.distance(to, j) +
                    p.lowerbound(to, pred);
                QNode* current = pQ.find(to);
                if(!current || newDistance < current->key)
                {
                    pQ.changeKey(QNode(newDistance, to), to);
                    pred.insert(to, j);
                }
            }
        }
    }
};

template<typename PROBLEM, typename STATE_ID = unsigned long long>
struct RecursiveBestFirstSearch
{
    Stack<STATE_ID> pred;
    PROBLEM& p;
    enum{SUCCESS = -1};
    typedef KVPair<double, STATE_ID> INFO;
    double work(INFO state, double alternative, double pathCost)
    {
        if(p.isGoal(state.value, pred)) return SUCCESS;
        Vector<STATE_ID> next = p.nextStates(state.value, pred);
        if(next.getSize() == 0) return numeric_limits<double>::max();
        Heap<INFO, KVComparator<double, STATE_ID> > children;
        for(int i = 0; i < next.getSize(); ++i)
            children.insert(INFO(max(state.key, pathCost +
                p.distance(state.value, next[i]) +
                p.lowerbound(next[i], pred)), next[i]));
        for(;;)
        {
            INFO best = children.deleteMin();
            if(best.key > alternative) return best.key;
            pred.push(best.value);
            best.key = work(best, children.isEmpty() ?
                alternative : min(children.getMin().key, alternative),
                pathCost + p.distance(state.value, best.value));
            if(best.key == SUCCESS) return SUCCESS;
            children.insert(best);
            pred.pop();
        }
    }
    RecursiveBestFirstSearch(PROBLEM& theProblem): p(theProblem)
    {
        pred.push(p.start());
        work(KVPair<double, int>(0.0, p.start()),
            numeric_limits<double>::max(), 0);
    }
};

template<typename PROBLEM> void branchAndBound(PROBLEM& p)
{
    if(!p.processSolution())
    {
        Vector<KVPair<double, typename PROBLEM::Move> > moves =
            p.generateMoves();
        quickSort(moves.getArray(), 0, moves.getSize() - 1,
            KVComparator<double, typename PROBLEM::Move>());
        for(int i = 0; i < moves.getSize(); ++i)
        {
            p.move(moves[i].value);
            branchAndBound(p);
            p.undoMove(moves[i].value);
        }
    }
}

}
#endif
