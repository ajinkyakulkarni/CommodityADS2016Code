#ifndef GRAPH_H
#define GRAPH_H
#include "../Utils/Vector.h"
#include "../Utils/Debug.h"
#include "../RandomNumberGeneration/Random.h"
#include "../Utils/Queue.h"
#include "../Utils/UnionFind.h"
#include <cassert>
#include "../Heaps/Heap.h"
using namespace std;
namespace igmdk{

template<typename EDGE_DATA> class GraphAA
{
    struct Edge
    {
        int to;
        EDGE_DATA edgeData;
        Edge(int theTo, EDGE_DATA const& theEdgeData): to(theTo),
            edgeData(theEdgeData) {}
    };
    Vector<Vector<Edge> > vertices;
public:
    GraphAA(int initialSize = 0): vertices(initialSize) {}
    int nVertices()const{return vertices.getSize();}
    int nEdges(int v)const{return vertices[v].getSize();}
    class AdjacencyIterator
    {
        Vector<Edge> const* edges;
        int j;
    public:
        AdjacencyIterator(GraphAA const& g, int v, int theJ):
            edges(&g.vertices[v]), j(theJ){}
        AdjacencyIterator& operator++()
        {
            ++j;
            return *this;
        }
        int to(){return (*edges)[j].to;}
        EDGE_DATA const& data(){return (*edges)[j].edgeData;}
        bool operator!=(AdjacencyIterator const& rhs){return j != rhs.j;}
    };
    AdjacencyIterator begin(int v)const
        {return AdjacencyIterator(*this, v, 0);}
    AdjacencyIterator end(int v)const
        {return AdjacencyIterator(*this, v, nEdges(v));}
    void addVertex(){vertices.append(Vector<Edge>());}
    void addEdge(int from, int to, EDGE_DATA const& edgeData = EDGE_DATA())
    {
        assert(to >= 0 && to < vertices.getSize());
        vertices[from].append(Edge(to, edgeData));
    }
    void addUndirectedEdge(int from, int to,
        EDGE_DATA const& edgeData = EDGE_DATA())
    {
        addEdge(from, to, edgeData);
        addEdge(to, from, edgeData);
    }
};

template<typename GRAPH> GRAPH reverse(GRAPH const& g)
{
    GRAPH result(g.nVertices());
    for(int i = 0; i < g.nVertices(); ++i)
        for(typename GRAPH::AdjacencyIterator j = g.begin(i);
            j != g.end(i); ++j) result.addEdge(j.to(), i, j.data());
    return result;
}

template<typename GRAPH>
GRAPH randomDirectedGraph(int vertices, int edgesPerVertex)
{
    assert(edgesPerVertex <= vertices);
    GRAPH g(vertices);
    for(int i = 0; i < vertices; ++i)
    {
        Vector<int> edges = GlobalRNG.sortedSample(edgesPerVertex, vertices);
        for(int j = 0; j < edgesPerVertex; ++j) g.addEdge(i, edges[i]);
    }
    return g;
}

template<typename GRAPH, typename ACTION> void DFSComponent(GRAPH const& g,
    int source, Vector<bool>& visited, ACTION& a = ACTION())
{
    Stack<int> vertexStack;
    typedef typename GRAPH::AdjacencyIterator ITER;
    Stack<ITER> nextStack;
    vertexStack.push(source);
    nextStack.push(g.begin(source));
    while(!vertexStack.isEmpty())
    {
        int v = vertexStack.getTop();
        ITER& j = nextStack.getTop();
        if(j != g.end(v))
        {
            if(visited[j.to()]) a.nonTreeEdge(j.to());
            else
            {
                a.treeEdge(j.to());
                visited[j.to()] = true;
                vertexStack.push(j.to());
                nextStack.push(g.begin(j.to()));
            }
            ++j;
        }
        else
        {
            vertexStack.pop();
            nextStack.pop();
            if(!vertexStack.isEmpty())
                a.backwardEdge(vertexStack.getTop());
        }
    }
}

template<typename GRAPH, typename ACTION> void DFS(GRAPH const& g,
    ACTION& a = ACTION())
{
    Vector<bool> visited(g.nVertices(), false);
    for(int i = 0; i < g.nVertices(); ++i) if(!visited[i])
    {
        a.source(i);
        visited[i] = true;
        DFSComponent(g, i, visited, a);
    }
}

struct ConnectedComponentAction
{
    Vector<Vector<int> > components;
    void source(int v)
    {
        components.append(Vector<int>());
        treeEdge(v);
    }
    void treeEdge(int v){components.lastItem().append(v);}
    void nonTreeEdge(int v){}
    void backwardEdge(int v){}
};
template<typename GRAPH>
Vector<Vector<int> > connectedComponents(GRAPH const& g)
{
    ConnectedComponentAction a;
    DFS(g, a);
    return a.components;
}

struct TopologicalSortAction
{
    int k, last;
    Vector<int> ranks;
    bool hasCycle;
    TopologicalSortAction(int nVertices): k(nVertices), last(-1),
        ranks(k, -1), hasCycle(false) {}
    void source(int v){treeEdge(v);}
    void treeEdge(int v){last = v;}
    void nonTreeEdge(int v){if(ranks[v] == -1) hasCycle = true;}
    void backwardEdge(int v)
    {
        if(last != -1)
        {
            ranks[last] = --k;
            last = -1;
        }
        ranks[v] = --k;
    }
};
template<typename GRAPH> Vector<int> topologicalSort(GRAPH const& g)
{//empty result means presence of cycle
    TopologicalSortAction a(g.nVertices());
    DFS(g, a);
    if(a.hasCycle) a.ranks = Vector<int>();
    return a.ranks;
}

template<typename GRAPH> Vector<int> BFS(GRAPH& g, int source)
{
    Vector<int> distances(g.nVertices(), -1);
    Queue<int> q(g.nVertices());
    distances[source] = 0;
    q.push(source);
    while(!q.isEmpty())
    {
        int i = q.pop();
        for(typename GRAPH::AdjacencyIterator j = g.begin(i); j != g.end(i);
            ++j) if(distances[j.to()] == -1)
            {
                distances[j.to()] = distances[i] + 1;
                q.push(j.to());
            }
    }
    return distances;
}

struct StaticGraph
{
    //does not need edge and vertex data, easy to iterate
    //inherently directed, for indirection needs coordination
    //by caller
    Vector<int> edges, starts;
    StaticGraph(){starts.append(0);}
    void addVertex()
    {
        starts.append(starts.lastItem());
        starts[starts.getSize()-2] = 0;
    }
    void addEdgeToLastVertex(int to)
    {
        edges.append(to);
        ++starts.lastItem();
    }
    int getSize(){return starts.getSize()-1;}
};

template<typename GRAPH> Vector<int> MST(GRAPH& g)
{
    Vector<int> mst(g.nVertices(), -1);
    typedef KVPair<double, int> QNode;
    IndexedArrayHeap<QNode, KVComparator<double, int> > pQ;
    for(int i = 0; i < g.nVertices(); ++i)
        pQ.insert(QNode(numeric_limits<double>::max(), i), i);
    while(!pQ.isEmpty())
    {
        int i = pQ.deleteMin().value;
        for(typename GRAPH::AdjacencyIterator j = g.begin(i); j != g.end(i);
            ++j)
        {
            QNode* current = pQ.find(j.to());
            if(current && j.data() < current->key)
            {
                pQ.changeKey(QNode(j.data(), j.to()), j.to());
                mst[j.to()] = i;
            }
        }
    }
    return mst;
}

template<typename GRAPH>
Vector<int> ShortestPath(GRAPH& g, int from, int dest = -1)
{
    assert(from >= 0 && from < g.nVertices());
    Vector<int> pred(g.nVertices(), -1);
    typedef KVPair<double, int> QNode;
    IndexedArrayHeap<QNode, KVComparator<double, int> > pQ;
    for(int i = 0; i < g.nVertices(); ++i) pQ.insert(
        QNode(i == from ? 0 : numeric_limits<double>::max(), i), i);
    while(!pQ.isEmpty() && pQ.getMin().value != dest)
    {
        int i = pQ.getMin().value;
        double dj = pQ.deleteMin().key;
        for(typename GRAPH::AdjacencyIterator j = g.begin(i); j != g.end(i);
            ++j)
        {
            double newDistance = dj + j.data();
            QNode* current = pQ.find(j.to());
            if(current && newDistance < current->key)
            {
                pQ.changeKey(QNode(newDistance, j.to()), j.to());
                pred[j.to()] = i;
            }
        }
    }
    return pred;
}

template<typename GRAPH> struct BellmanFord
{
    int v;
    Vector<double> distances;
    Vector<int> pred;
    bool hasNegativeCycle;
    bool findNegativeCycle()
    {
        UnionFind uf(v);
        for(int i = 0; i < v; ++i)
        {
            int p = pred[i];
            if(p != -1)
            {
                if(uf.areEquivalent(i, p)) return true;
                uf.join(i, p);
            }
        }
        return false;
    }
    BellmanFord(GRAPH& g, int from): v(g.nVertices()), pred(v, -1),
        distances(v, numeric_limits<double>::max()), hasNegativeCycle(false)
    {
        assert(from >= 0 && from < v);
        Vector<bool> onQ(v, 0);
        Queue<int> queue;
        int cost = 0;
        distances[from] = 0;
        onQ[from] = true;
        queue.push(from);
        do
        {
            int i = queue.pop();
            onQ[i] = false;
            for(typename GRAPH::AdjacencyIterator j = g.begin(i);
                j != g.end(i); ++j)
            {
                double newDistance = distances[i] + j.data();
                if(newDistance < distances[j.to()])
                {
                    distances[j.to()] = newDistance;
                    pred[j.to()] = i;
                    if(!onQ[j.to()])
                    {
                        queue.push(j.to());
                        onQ[j.to()] = true;
                    }
                }
                if(++cost % v == 0) hasNegativeCycle = findNegativeCycle();
            }
        }while(!queue.isEmpty() && !hasNegativeCycle);
    }
};

struct FlowData
{
    int from;
    double flow, capacity, cost;//cost used only for min flow
    FlowData(int theFrom, double theCapacity, double theCost = 0):
        from(theFrom), capacity(theCapacity), flow(0), cost(theCost) {}
    double capacityTo(int v){return v == from ? flow : capacity - flow;}
    void addFlowTo(int v, double change)
        {flow += change * (v == from ? -1 : 1);}
};
template<typename GRAPH> class ShortestAugmentingPath
{
    int v;
    Vector<int> path, pred;
    bool hasAugmentingPath(GRAPH& g, Vector<FlowData>& fedges, int from,
        int to)
    {
        for(int i = 0; i < v; ++i) pred[i] = -1;
        Queue<int> queue;
        queue.push(pred[from] = from);
        while(!queue.isEmpty())
        {
            int i = queue.pop();
            for(typename GRAPH::AdjacencyIterator j = g.begin(i);
                j != g.end(i); ++j)
                if(pred[j.to()] == -1 &&
                    fedges[j.data()].capacityTo(j.to()) > 0)
                {
                    path[j.to()] = j.data();
                    pred[j.to()] = i;
                    queue.push(j.to());
                }
        }
        return pred[to] != -1;
    }
    bool hasMinCostAugmentingPath(GRAPH& g, Vector<FlowData>& fedges,
        int from, int to)
    {
        if(total >= neededFlow) return false;
        GraphAA<double> costGraph(v);
        for(int i = 0; i < v; ++i)
            for(typename GRAPH::AdjacencyIterator j = g.begin(i);
                j != g.end(i); ++j)
                if(fedges[j.data()].capacityTo(j.to()) > 0)
                    costGraph.addEdge(i, j.to(), fedges[j.data()].cost);
        pred = ShortestPath(costGraph, from, to);
        for(int i = to; pred[i] != -1; i = pred[i])
            for(typename GRAPH::AdjacencyIterator j = g.begin(pred[i]);
                j != g.end(pred[i]); ++j)
                if(j.to() == i) path[j.to()] = j.data();
        return pred[to] != -1;
    }
public:
    double total, neededFlow;
    ShortestAugmentingPath(GRAPH& g, Vector<FlowData>& fedges, int from,
        int to, double theNeededFlow = 0): v(g.nVertices()), total(0),
        path(v, -1), pred(v, -1), neededFlow(theNeededFlow)
    {
        assert(from >= 0 && from < v && to >= 0 && to < v);
        while(neededFlow > 0 ? hasMinCostAugmentingPath(g, fedges, from, to)
            : hasAugmentingPath(g, fedges, from, to))
        {
            double increment = numeric_limits<double>::max();
            for(int j = to; j != from; j = pred[j])
                increment = min(increment, fedges[path[j]].capacityTo(j));
            if(neededFlow > 0)increment = min(increment, neededFlow - total);
            for(int j = to; j != from; j = pred[j])
                fedges[path[j]].addFlowTo(j, increment);
            total += increment;
        }
    }
};

Vector<pair<int, int> > bipartiteMatching(int n, int m,
    Vector<pair<int, int> > const& allowedMatches)
{//v = n + m + 2, e = n + m + allowedMatches.getSize(), time is O(ve)
    GraphAA<int> sp(n + m + 2);
    Vector<FlowData> data;
    for(int i = 0; i < allowedMatches.getSize(); ++i)
    {
        data.append(FlowData(allowedMatches[i].first, 1));
        sp.addUndirectedEdge(allowedMatches[i].first,
            allowedMatches[i].second, i);
    }
    int source = n + m, sink = source + 1;
    for(int i = 0; i < source; ++i)
    {
        int from = i, to = sink;
        if(i < n)
        {
            from = source;
            to = i;
        }
        data.append(FlowData(from, 1));
        sp.addUndirectedEdge(from, to, i + allowedMatches.getSize());
    }
    ShortestAugmentingPath<GraphAA<int> > dk(sp, data, source, sink);
    Vector<pair<int, int> > result;
    for(int i = 0; i < allowedMatches.getSize(); ++i)
        if(data[i].flow > 0) result.append(allowedMatches[i]);
    return result;
}

Vector<pair<int, int> > assignmentProblem(int n, int m,
    Vector<pair<pair<int, int>, double> > const& allowedMatches)
{//v = n + m + 2, e = n + m + allowedMatches.getSize()
    GraphAA<int> sp(n + m + 2);
    Vector<FlowData> data;
    for(int i = 0; i < allowedMatches.getSize(); ++i)
    {
        data.append(FlowData(allowedMatches[i].first.first, 1,
            allowedMatches[i].second));
        sp.addUndirectedEdge(allowedMatches[i].first.first,
            allowedMatches[i].first.second, i);
    }
    int source = n + m, sink = source + 1;
    for(int i = 0; i < source; ++i)
    {
        int from = i, to = sink;
        if(i < n)
        {
            from = source;
            to = i;
        }
        data.append(FlowData(from, 1, 0));
        sp.addUndirectedEdge(from, to, i + allowedMatches.getSize());
    }
    ShortestAugmentingPath<GraphAA<int> > dummy(sp, data, source, sink,
        min(n, m));
    Vector<pair<int, int> > result;
    for(int i = 0; i < allowedMatches.getSize(); ++i)
        if(data[i].flow > 0) result.append(allowedMatches[i].first);
    return result;
}

Vector<int> stableMatching(Vector<Vector<int> > const& womenOrders,
    Vector<Vector<int> > const& menScores)
{
    int n = womenOrders.getSize(), m = menScores.getSize();
    assert(n <= m);
    Stack<int> unassigned;//any list type will do
    for(int i = 0; i < n; ++i) unassigned.push(i);
    Vector<int> currentMan(m, -1), nextWoman(n, 0);
    while(!unassigned.isEmpty())
    {
        int man = unassigned.pop(), woman, current;
        do//man finds best woman that prefers him to her current partner
        {
            woman = nextWoman[man]++;
            current = currentMan[woman];
        }while(current != -1 &&
            menScores[woman][man] <= menScores[woman][current]);
        if(current != -1) unassigned.push(current);
        currentMan[woman] = man;
    }
    return currentMan;
}

}
#endif
