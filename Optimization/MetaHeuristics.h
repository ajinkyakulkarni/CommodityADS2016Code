#ifndef METAHEURISTICS_H
#define METAHEURISTICS_H
#include "../RandomNumberGeneration/Random.h"
namespace igmdk{

template<typename SOLUTION> void localSearch(SOLUTION& s, long long maxMoves,
    int maxStall)
{
    for(int i = 0; maxMoves-- && i < maxStall; ++i)
        if(s.proposeMove() >= 0)
        {
            i = -1;
            s.applyMove();
        }
}

template<typename SOLUTION> void simulatedAnnealing(SOLUTION& s, double T,
    double coolingFactor, long long maxMoves)
{
    while(maxMoves--)
    {
        double change = s.proposeMove();
        if(change > -T * GlobalRNG.exponential(1)) s.applyMove();
        T *= coolingFactor;
    }
}

template<typename SOLUTION>
void iteratedLocalSearch(SOLUTION& s, long long maxBigMoves)
{
    while(maxBigMoves--)
    {
        s.localSearchBest();
        s.updateBest();
        s.bigMove();
    }
}

}
#endif
