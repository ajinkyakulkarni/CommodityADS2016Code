#ifndef STRING_ALGS_H
#define STRING_ALGS_H
#include "../Utils/Utils.h"
#include "../Utils/Vector.h"
#include "../Utils/GCFreelist.h"
#include "../Utils/Stack.h"
#include "../Graphs/Graph.h"
namespace igmdk{

enum{ALPHABET_SIZE = 1 << numeric_limits<unsigned char>::digits};

template<typename VECTOR> bool matchesAt(int position, VECTOR text,
    VECTOR pattern, int patternSize)
{
    int i = 0;
    while(i < patternSize && pattern[i] == text[i + position]) ++i;
    return i == patternSize;
}

class Horspool
{
    int shift[ALPHABET_SIZE], position, textSize, patternSize;
    unsigned char *text, *pattern;
public:
    Horspool(unsigned char* theText, int theTextSize,
        unsigned char* thePattern, int thePatternSize): position(0),
        textSize(theTextSize), text(theText), patternSize(thePatternSize),
        pattern(thePattern)
    {
        for(int i = 0; i < ALPHABET_SIZE; ++i) shift[i] = patternSize;
        for(int i = 0; i < patternSize - 1; ++i)
            shift[pattern[i]] = patternSize - 1 - i;
    }
    int findNext()
    {
        while(position + patternSize <= textSize)
        {
            int result = position;
            position += shift[text[position + patternSize - 1]];
            if(matchesAt(result, text, pattern, patternSize)) return result;
        }
        return -1;
    }
};

struct Q3Hash
{
    int size;
    int operator()(unsigned char* substring, int position)
    {//assumes q = 3;
        substring += position;
        unsigned char result = substring[0];
        result += substring[1] << 1;
        result += substring[2] << 2;
        return result % size;
    }
};
template<typename VECTOR, typename HASH> class HashQ
{
    int position, textSize, patternSize, q;
    Vector<int> shift;
    VECTOR const &text, &pattern;
    HASH h;
public:
    HashQ(VECTOR const& theText, int theTextSize, VECTOR const& thePattern,
        int thePatternSize, int theQ = 3, HASH theHash = HASH()):
        position(0), textSize(theTextSize), text(theText), q(theQ),
        pattern(thePattern), patternSize(thePatternSize),
        shift(patternSize/q), h(theHash)
    {
        h.size = shift.getSize();
        assert(patternSize >= q);
        int temp = patternSize - q;
        for(int i = 0; i < shift.getSize(); ++i) shift[i] = temp + 1;
        for(int i = 0; i < temp; ++i) shift[h(pattern, i)] = temp - i;
    }
    int findNext()
    {
        while(position + patternSize <= textSize)
        {
            int result = position;
            position += shift[h(text, position + patternSize - q)];
            if(matchesAt(result, text, pattern, patternSize))return result;
        }
        return -1;
    }
};

template<typename HASH = Q3Hash> class WuManber
{
    unsigned char* text;
    Vector<pair<unsigned char*, int> >const& patterns;
    int shift[ALPHABET_SIZE], position, q, minPatternSize, textSize;
    Vector<int> candidates[ALPHABET_SIZE];
    HASH h;
    void findMatches(int h, Vector<int>& matches)
    {
        for(int i = 0; i < candidates[h].getSize(); ++i)
        {
            int j = candidates[h][i], patternISize = patterns[j].second;
            if(position + patternISize <= textSize && matchesAt(position,
                text, patterns[j].first, patternISize))  matches.append(j);
        }
    }
public:
    WuManber(unsigned char* theText, int theTextSize,
        Vector<pair<unsigned char*, int> >const& thePatterns,
        HASH theHash = HASH()): text(theText),
        patterns(thePatterns), position(0), h(theHash), q(3),
        textSize(theTextSize), minPatternSize(numeric_limits<int>::max())
    {
        h.size = ALPHABET_SIZE;
        for(int i = 0; i < patterns.getSize(); ++i)
            minPatternSize = min(patterns[i].second, minPatternSize);
        assert(minPatternSize >= q);
        int temp = minPatternSize - q;
        for(int i = 0; i < ALPHABET_SIZE; ++i) shift[i] = temp + 1;
        for(int j = 0; j < patterns.getSize(); ++j)
            for(int i = 0; i < temp + 1; ++i)
            {
                int hi = h(patterns[j].first, i);
                if(i == temp) candidates[hi].append(j);
                else shift[hi] = min(temp - i, shift[hi]);
            }
    }
    int findNext(Vector<int>& matches)
    {
        while(position + minPatternSize <= textSize)
        {
            int hi = h(text, position + minPatternSize - q),
                result = position;
            findMatches(hi, matches);
            position += shift[hi];
            if(matches.getSize() > 0) return result;
        }
        return -1;
    }
};

struct RegularExpressionMatcher
{
    string re;
    int m;
    GraphAA<bool> g;
    RegularExpressionMatcher(string const& theRe): re(theRe), m(re.length()),
        g(m + 1)
    {
        Stack<int> clauses;
        for(int i = 0; i < m; ++i)
        {
            int clauseStart = i;
            if(re[i] == '(' || re[i] == '|') clauses.push(i);
            else if(re[i] == ')')
            {
                int clauseOp = clauses.pop();
                if(re[clauseOp] == '|')
                {
                    clauseStart = clauses.pop();
                    g.addEdge(clauseStart, clauseOp+1);
                    g.addEdge(clauseOp, i);
                }
                else clauseStart = clauseOp;
            }
            if(i < m - 1 && re[i + 1]=='*')
                g.addUndirectedEdge(clauseStart, i + 1);
            if(re[i] == '(' || re[i] == '*' || re[i] == ')')
                g.addEdge(i, i + 1);
        }
    }

    struct VisitedAction
    {
        Vector<bool> visited;
        VisitedAction(int nVertices): visited(nVertices, false) {}
        void treeEdge(int v){}
        void nonTreeEdge(int v){}
        void backwardEdge(int v){}
    };
    Vector<int> findActiveStates(Vector<int> sources = Vector<int>())
    {
        VisitedAction a(g.nVertices());
        if(sources.getSize() == 0) DFSComponent(g, 0, a.visited, a);
        else for(int i = 0; i < sources.getSize(); ++i)
                if(!a.visited[sources[i]])
                {
                    a.visited[sources[i]] = true;
                    DFSComponent(g, sources[i], a.visited, a);
                }
        Vector<int> activeStates;
        for(int i = 0; i < a.visited.getSize(); ++i)
            if(a.visited[i]) activeStates.append(i);
        return activeStates;
    }

    bool matches(string const& text)
    {
        Vector<int> activeStates = findActiveStates();
        for(int i = 0; i < text.length(); ++i)
        {
            Vector<int> match;
            for(int j = 0; j < activeStates.getSize(); ++j)
                if(activeStates[j] < m && re[activeStates[j]] == text[i])
                    match.append(activeStates[j] + 1);
            activeStates = findActiveStates(match);
        }
        for(int j = 0; j < activeStates.getSize(); ++j)
            if(activeStates[j] == m) return true;
        return false;
    }
};

class ShiftAnd
{
    unsigned char *pattern, *text;
    int textSize, patternSize, position;
    unsigned long long charPos[ALPHABET_SIZE], state;
public:
    ShiftAnd(unsigned char* theText, int theTextSize,
        unsigned char* thePattern, int thePatternSize): position(0),
        textSize(theTextSize), text(theText), patternSize(thePatternSize),
        pattern(thePattern), state(0)
    {
        assert(patternSize <= numeric_limits<unsigned long long>::digits);
        for(int i = 0; i < ALPHABET_SIZE; ++i) charPos[i] = 0;
        for(int i = 0; i < patternSize; ++i)
            charPos[pattern[i]] |= 1ull << i;
    }
    int findNext()
    {
        unsigned long long limit = 1ull << (patternSize - 1);
        while(position < textSize)
        {
            state = ((state << 1) | 1) & charPos[text[position++]];
            if(state & limit) return position - patternSize;
        }
        return -1;
    }
};

class ShiftAndExtended
{
    unsigned char *pattern, *text;
    int textSize, patternSize, position;
    unsigned long long *charPos, P, L, S, state;
public:
    ShiftAndExtended(unsigned char* theText, int theTextSize,
        int thePatternSize, unsigned long long* theCharPos,
        unsigned long long optionalPos, unsigned long long repeatedPos):
        position(0), textSize(theTextSize), text(theText),
        patternSize(thePatternSize), S(repeatedPos), state(0)
    {
        assert(patternSize <= numeric_limits<unsigned long long>::digits);
        unsigned long long sides = optionalPos ^ (optionalPos >> 1);
        P = (optionalPos >> 1) & sides;
        L = optionalPos & sides;
    }
    int findNext()
    {
        unsigned long long limit = 1ull << (patternSize - 1);
        while(position < textSize)
        {
            state = (((state << 1) | 1) | (state & S)) &
                charPos[text[position++]];
            state |= L ^ ~(L - (state & P));
            if(state & limit) return position - patternSize;
        }
        return -1;
    }
};

struct Edit
{
    Edit* prev;
    int position;
    bool isInsert;
};
template<typename CHAR> void extendDiagonal(int diagonal,
    Vector<int>& frontierX, Vector<Edit*>& edits, Vector<CHAR>const& a,
    Vector<CHAR>const& b, Freelist<Edit>& freelist)
{
    int x = max(frontierX[diagonal - 1] + 1, frontierX[diagonal + 1]),
        y = x - (diagonal - 1 - a.getSize());
    if(x != -1 || y != -1)
    {
        bool isInsert = x != frontierX[diagonal + 1];
        edits[diagonal] = new(freelist.allocate())Edit();
        edits[diagonal]->isInsert = isInsert;
        edits[diagonal]->prev = edits[diagonal + (isInsert ? -1 : 1)];
        edits[diagonal]->position = isInsert ? x : y;
    }
    while(y + 1 < a.getSize() && x + 1 < b.getSize() && a[y + 1] == b[x + 1])
    {
        ++y;
        ++x;
    }
    frontierX[diagonal] = x;
}
template<typename CHAR>
Vector<Edit> DiffInternal(Vector<CHAR>const& a, Vector<CHAR>const& b)
{
    int M = a.getSize(), N = b.getSize(), size = M + N + 3, d = N + 1;
    assert(M <= N);//a must be shorter then b
    Vector<int> frontierX(size, -2);
    Vector<Edit*> edits(size, 0);
    Freelist<Edit> freelist;
    int p = 0;
    for(; frontierX[d] < N - 1; ++p)
    {
        for(int diagonal = M + 1 - p; diagonal < d; ++diagonal)
            extendDiagonal<CHAR>(diagonal, frontierX, edits, a, b, freelist);
        for(int diagonal = d + p; diagonal >= d; --diagonal)
            extendDiagonal<CHAR>(diagonal, frontierX, edits, a, b, freelist);
    }
    Vector<Edit> result;
    for(Edit* link = edits[d]; link; link = link->prev) result.append(*link);
    result.reverse();
    //distance is N - M + 2 * (p - 1) with p - 1 deletions
    assert(result.getSize() == N - M + 2 * (p - 1));
    return result;
}
template<typename CHAR>
Vector<Edit> Diff(Vector<CHAR>const& a, Vector<CHAR>const& b)
{//edits needed to get a into b
    bool change = a.getSize() > b.getSize();
    Vector<Edit> result = change ? DiffInternal(b, a) : DiffInternal(a, b);
    for(int i = 0, offset = 0; i < result.getSize(); ++i)
    {
        if(change) result[i].isInsert = !result[i].isInsert;
        if(result[i].isInsert) ++offset;
        else result[i].position += offset--;
    }
    return result;
}

}//end namespace
#endif

