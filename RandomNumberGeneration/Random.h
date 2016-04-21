#ifndef RANDOM_H
#define RANDOM_H
#include <ctime>
#include <cassert>
#include <cmath>
#include <limits>
#include "../Utils/Utils.h"
#include "../Utils/Vector.h"
#include "../Utils/Stack.h"
using namespace std;
namespace igmdk{

class Xorshift
{
    unsigned int state;
    enum{PASSWORD = 19870804};
public:
    Xorshift(unsigned int seed = time(0) ^ PASSWORD)
    {
        assert(numeric_limits<unsigned int>::digits == 32);
        state = seed ? seed : PASSWORD;
    }
    static unsigned int transform(unsigned int x)
    {
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        return x;
    }
    unsigned int next(){return state = transform(state);}
    double uniform01(){return 2.32830643653869629E-10 * next();}
};

class ImprovedXorshift
{
    unsigned int state1, state2;
    enum{PASSWORD = 19870804};
public:
    ImprovedXorshift(unsigned int seed = time(0) ^ PASSWORD)
    {
        assert(numeric_limits<unsigned int>::digits == 32);
        state1 = seed ? seed : PASSWORD;
        state2 = state1;
    }
    unsigned int next()
    {
        state1 ^= state1 << 13;
        state1 ^= state1 >> 17;
        state1 ^= state1 << 5;
        state2 = state2 * 69069U + 1234567U;
        return state1 + state2;
    }//may return 0
    double uniform01(){return 2.32830643653869629E-10 * max(1u, next());}
};

class QualityXorshift64
{
    unsigned long long state;
    enum{PASSWORD = 19870804};
public:
    QualityXorshift64(unsigned long long seed = time(0) ^ PASSWORD)
    {
        assert(numeric_limits<unsigned long long>::digits == 64);
        state = seed ? seed : PASSWORD;
    }
    static unsigned long long transform(unsigned long long x)
    {
        x ^= x << 21;
        x ^= x >> 35;
        x ^= x << 4;
        return x * 2685821657736338717ull;
    }
    unsigned long long next(){return state = transform(state);}
    double uniform01(){return 5.42101086242752217E-20 * next();}
};

struct MRG32k3a
{
    enum{PASSWORD = 19870804};
    static long long const m1 = 4294967087ll, m2 = 4294944443ll;
    long long s10, s11, s12, s20, s21, s22;
    void reduceAndUpdate(long long c1, long long c2)
    {
        if(c1 < 0) c1 = m1 - (-c1 % m1);
        else c1 %= m1;
        if(c2 < 0) c2 = m2 - (-c2 % m2);
        else c2 %= m2;
        s10 = s11; s11 = s12; s12 = c1;
        s20 = s21; s21 = s22; s22 = c2;
    }
public:
    unsigned int next()
    {
        long long c1 = (1403580 * s11 - 810728 * s10),
            c2 = (527612 * s22 - 1370589 * s20);
        reduceAndUpdate(c1, c2);
        return (c1 <= c2 ? m1 : 0) + c1 - c2;
    }
    //s1(0-2) and s2(0-2) must be respectively < m1 and m2 and not all 0
    MRG32k3a(): s10(max(time(0) ^ PASSWORD, 1l) % m2), s11(0), s12(0),
        s20(s10), s21(0), s22(0){}
    double uniform01(){return next()/(m1 + 1.0);}
    void jumpAhead()
    {
        const long long A1p76[3][3] = {
            {  82758667u, 1871391091u, 4127413238u},//for s10
            {3672831523u,   69195019u, 1871391091u},//for s11
            {3672091415u, 3528743235u,   69195019u}},//for s12
            A2p76[3][3] = {
                {1511326704u, 3759209742u, 1610795712u},//for s20
                {4292754251u, 1511326704u, 3889917532u},//for s21
                {3859662829u, 4292754251u, 3708466080u}};//for s22
        long long s1[3] = {s10, s11, s12}, s2[3] = {s20, s21, s22};
        for(int i = 0; i < 3; ++i)
        {
            long long c1 = 0, c2 = 0;
            for(int j = 0; j < 3; ++j)
            {
                c1 += s1[j] * A1p76[i][j];
                c2 += s2[j] * A2p76[i][j];
            }
            reduceAndUpdate(c1, c2);
        }
    }
};

struct ARC4
{
    unsigned char sBox[256], i, j;
    enum{PASSWORD = 19870804};
    void construct(unsigned char* seed, int length)
    {
        j = 0;
        for(int k = 0; k < 256; ++k) sBox[k] = k;
        for(int k = 0; k < 256; ++k)
        {//different from the random permutation algorithm
            j += sBox[k] + seed[k % length];
            swap(sBox[k], sBox[j]);
        }
        i = j = 0;
        for(int dropN = 1024; dropN > 0; dropN--) nextByte();
    }
    ARC4(unsigned long long seed = time(0) ^ PASSWORD)
        {construct((unsigned char*)&seed, sizeof(seed));}
    //for crytographic initialization from a long seed
    ARC4(unsigned char* seed, int length){construct(seed, length);}
    unsigned char nextByte()
    {
        j += sBox[++i];
        swap(sBox[i], sBox[j]);
        return sBox[(sBox[i] + sBox[j]) % 256];
    }
    unsigned long long next()
    {
        unsigned long long result = 0;
        for(int k = 0; k < sizeof(result); ++k)
            result |= ((unsigned long long)nextByte()) << (8 * k);
        return result;
    }
    double uniform01(){return 5.42101086242752217E-20 * max(1ull, next());}
};

template<typename GENERATOR = QualityXorshift64> struct Random
{
    GENERATOR random;
    enum{PASSWORD = 19870804};
    Random(unsigned long long seed = time(0) ^ PASSWORD): random(seed){}
    unsigned long long next(){return random.next();}
    unsigned long long mod(unsigned long long n)
    {
        assert(n > 0);
        return next() % n;
    }
    long long inRange(long long a, long long b){return a + mod(b - a + 1);}
    double uniform01(){return random.uniform01();}
    bool bernoulli(double p){return uniform01() <= p;}
    int binomial(double p, int n)
    {
        int result = 0;
        for(int i = 0; i < n; ++i) result += bernoulli(p);
        return result;
    }
    int geometric(double p)
    {
        assert(p > 0);
        int result = 0;
        while(!bernoulli(p)) ++result;
        return result;
    }
    int geometric05(){return rightmost0Count(next());}
    int poisson(double l)
    {
        assert(l > 0);
        int result = -1;
        for(double p = 1; p > exp(-l); p *= uniform01()) ++result;
        return result;
    }
    double uniform(double a, double b){return a + (b - a) * uniform01();}
    double exponential(double a){return -log(uniform01())/a;}
    static double PI(){return 3.1415926535897932384626433832795;}
    double cauchy(double m, double q)
        {return (tan((uniform01() - 0.5) * PI()) + m) * q;}
    double weibull1(double b){return pow(exponential(1), 1/b);}
    double normal01()
    {
        for(;;)
        {
            double a = 2 * uniform01() - 1, b = 2 * uniform01() - 1,
                c = a * a + b * b;
            if(c < 1)
            {
                double temp = sqrt(-2 * log(c)/c);
                return a * temp;//can return b * temp as 2nd iid sample
            }
        }
    }
    double normal(double m, double stdev){return m + stdev * normal01();}
    double logNormal(double m, double q){return exp(normal(m, q));}
    double gamma1(double b)
    {
        if(b >= 1)
        {
            for(double third = 1.0/3, d = b - third, x, v, u, xs;;)
            {
                do
                {
                    x = normal01();
                    v = 1 + x * third/sqrt(d);
                }while(v <= 0);
                v *= v * v; u = uniform01(), xs = x * x;
                if(u > 0.0331 * xs * xs || log(u) < xs/2 +
                    d * (1 - v + log(v))) return d * v;
            }
        }
        else
        {
            assert(b > 0);
            return pow(uniform01(), 1/b) * gamma1(b + 1);
        }
    }
    double erlang(double m, int k){return gamma1(k) * m/k;}
    double chiSquared(int k){return 2 * gamma1(k/2.0);}
    double t(int v){return sqrt(v/chiSquared(v)) * normal01();}
    double beta(double p, double q)
    {
        double G1 = gamma1(p);
        return G1/(G1 + gamma1(q));
    }
    double F(int v1, int v2)
        {return v2 * chiSquared(v1)/(v1 * chiSquared(v2));}
    double triangular01(double middle)
    {
        double u = uniform01();
        return sqrt(u <= middle ? middle * u : (1 - middle) * (1 - u));
    }
    double triangular(double a, double b, double c)
        {a + (b - a) * triangular01((c - a)/(b - a));}
    template<typename ITEM> void randomPermutation(ITEM* numbers, int size)
    {
        for(int i = 0; i < size; ++i)
            swap(numbers[i], numbers[inRange(i, size - 1)]);
    }
    Vector<int> sortedSample(int k, int n)
    {
        Vector<int> result;
        for(int considered = 0, selected = 0; selected < k; ++considered)
            if(bernoulli(double(k - selected)/(n - considered)))
            {
                result.append(considered);
                ++selected;
            }
        return result;
    }
    double uniformOrderStatistic(int i, int n){return beta(i, n - i + 1);}
};
Random<> GlobalRNG;

template<typename ITEM> void permuteDeterministically(ITEM* a, int size)
{
    Random<> drng(0);
    drng.randomPermutation(a, size);
}

class AliasMethod
{
    int n;
    Vector<int> aliases;
    Vector<double> wealth;
public:
    AliasMethod(Vector<double> const& probabilities):
        n(probabilities.getSize()), aliases(n, -1), wealth(n, 0)
    {
        Stack<int> smaller, greater;
        for(int i = 0; i < n; ++i)
        {//separate into poor and rich
            (wealth[i] = n * probabilities[i]) < 1 ?
                smaller.push(i) : greater.push(i);
        }
        while(!smaller.isEmpty() && !greater.isEmpty())
        {//reassign wealth until no poor remain
            int rich = greater.getTop(), poor = smaller.pop();
            aliases[poor] = rich;
            wealth[rich] -= 1 - wealth[poor];
            if(wealth[rich] < 1) smaller.push(greater.pop());
        }
    }
    int next()
    {//-1 check handles wealth round-off accumulation
        int x = GlobalRNG.mod(n);
        return GlobalRNG.uniform01() < wealth[x] || aliases[x] == -1 ?
            x : aliases[x];
    }
};

template<typename ITEM> class SumHeap
{
    Vector<ITEM> heap;
    int parent(int i){return (i - 1)/2;}
    int leftChild(int i){return 2 * i + 1;}
public:
    ITEM total(){return heap[0];}
    void add(ITEM value, int i = -1)
    {
        if(i == -1)
        {
            i = heap.getSize();
            heap.append(0);
        }
        for(;; i = parent(i))
        {
            heap[i] += value;
            if(i < 1) break;
        }
    }

    int find(ITEM value)
    {
        assert(0 <= value && value <= total());
        ITEM left = 0;
        for(int i = 0, c;; i = c)
        {
            c = leftChild(i);
            if(c >= heap.getSize()) return i;
            if(value > left + heap[c])
            {
                left += heap[c++];
                if(c >= heap.getSize() || value > left + heap[c]) return i;
            }
        }
    }
    int next(){return find(GlobalRNG.uniform01()*total());}

    ITEM cumulative(int i)
    {
        ITEM sum = heap[i];
        while(i > 0)
        {//add value of every left sibling if
            int last = i;
            i = parent(i);
            int l = leftChild(i);
            if(l != last) sum += heap[l];
        }
        return sum;
    }
    ITEM get(int i)
    {
        ITEM result = heap[i];
        int c = leftChild(i);
        if(c < heap.getSize())
        {
            result -= heap[c];
            if(++c < heap.getSize()) result -= heap[c];
        }
        return result;
    }
};

template<typename ITEM> struct ReservoirSampler
{
    int k, nProcessed;
    Vector<ITEM> selected;
    void processItem(ITEM const& item)
    {
        ++nProcessed;
        if(selected.getSize() < k) append(item);
        else
        {
            int kickedOut = GlobalRNG.mod(nProcessed);
            if(kickedOut < k) selected[kickedOut] = item;
        }
    }
    ReservoirSampler(int wantedSize): k(wantedSize), nProcessed(0){}
};

template<typename BIASED_COIN> class FairCoin
{
    BIASED_COIN biasedCoin;
    FairCoin(BIASED_COIN const& theBiasedCoin = BIASED_COIN()):
        biasedCoin(theBiasedCoin){}
    bool flip()
    {//HT is 1, TH is 0, ignore HH and TT
        bool flip1;
        do
        {
            flip1 = biasedCoin.flip();
        }while(flip1 == biasedCoin.flip());
        return flip1;
    }
};

void normalizeProbs(Vector<double>& probs)
{
    double sum = 0;
    for(int i = 0; i < probs.getSize(); ++i) sum += probs[i];
    probs *= 1/sum;
}

}
#endif
