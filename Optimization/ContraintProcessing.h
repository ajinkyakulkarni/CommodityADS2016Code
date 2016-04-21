#ifndef CONSTRAINT_PROCESSING_H
#define CONSTRAINT_PROCESSING_H
#include "../Utils/Queue.h"
#include "../Graphs/Graph.h"
#include "../Utils/Bitset.h"
#include "../HashTable/LinearProbingHashTable.h"
namespace igmdk{

template<typename CONSTRAINT> struct ConstraintGraph
{
    typedef GraphAA<CONSTRAINT> GRAPH;
    GRAPH g;
    Vector<Bitset<> > variables;
    void addVariable(int domain)
    {
        g.addVertex();
        variables.append(Bitset<>(domain));
    }
    void addConstraint(int v1, int v2, CONSTRAINT const& constraint)
    {
        assert(v1 != v2);
        g.addUndirectedEdge(v1, v2, constraint);
    }
    void disallow(int variable, int value)
        {variables[variable].set(value, false);}
    bool hasSolution(int variable){return !variables[variable].isZero();}

    bool isAllowed(int variable, int value, int otherVariable,
        CONSTRAINT const& constraint)
    {
        for(int i = 0; i < variables[otherVariable].getSize(); ++i)
            if(variables[otherVariable][i] && constraint.isAllowed(variable,
                value, otherVariable, i)) return true;
        return false;
    }
    bool revise(int variable1, int variable2, CONSTRAINT const& constraint)
    {
        bool changed = false;
        for(int i = 0; i < variables[variable1].getSize(); ++i)
            if(variables[variable1][i] &&
                !isAllowed(variable1, i, variable2, constraint))
            {
                disallow(variable1, i);
                changed = true;
            }
        return changed;
    }
    bool AC3Helper(int v, Queue<int>& q, Vector<bool>& onQ, bool isFirstPass)
    {
        onQ[v] = false;
        for(typename GRAPH::AdjacencyIterator i = g.begin(v);
            i != g.end(v); ++i)
        {
            int revisee = i.to(), against = v;
            if(isFirstPass) swap(revisee, against);
            if(revise(revisee, against, i.data()))
            {
                if(!hasSolution(revisee)) return false;
                if(!onQ[revisee]) q.push(revisee);
                onQ[revisee] = true;
            }
        }
        return true;
    }
    bool AC3()
    {
        Queue<int> q;
        Vector<bool> onQ(g.nVertices(), true);
        for(int j = 0; j < g.nVertices(); ++j)
            if(!AC3Helper(j, q, onQ, true)) return false;
        while(!q.isEmpty())if(!AC3Helper(q.pop(), q, onQ, false))
            return false;
        return true;
    }
};

struct AllDifferent
{
    LinearProbingHashTable<int, bool> variables;
    void addVariable(int variable){variables.insert(variable, true);}
    struct Handle
    {
        LinearProbingHashTable<int, bool>& variables;
        bool isAllowed(int variable, int value, int variable2, int value2)
            const
        {
            if(variables.find(variable) && variables.find(variable2))
                return value != value2;
            return true;
        }
        Handle(LinearProbingHashTable<int, bool>& theVariables):
            variables(theVariables) {}
    } handle;
    AllDifferent(): handle(variables) {}
};
}
#endif
