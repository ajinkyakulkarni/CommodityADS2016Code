
#ifndef STATISTICS_H
#define STATISTICS_H
#include "Random.h"
#include "../Utils/Sort.h"
#include "../Utils/Bits.h"
#include "../Heaps/Heap.h"
#include <cstdlib>
namespace igmdk{

double approxErf(double x)
{//for 0 <= x < inf, max error = 3e-7
    double a[6] = {0.0705230784, 0.0422820123, 0.0092705272, 0.0001520143,
        0.0002765672, 0.0000430638}, poly = 1, xPower = x;
    for(int i = 0; i < 6; ++i)
    {
        poly += a[i] * xPower;
        xPower *= x;
    }
    for(int i = 0; i < 4; ++i) poly *= poly;
    return 1 - 1/poly;
}
double approxNormalCDF(double x){return 0.5 + approxErf(x/sqrt(2))/2;}
double approxNormal2SidedConf(double x){return 2 * approxNormalCDF(x) - 1;}

struct NormalSummary
{
    double mean, variance;
    double stddev()const{return sqrt(variance);}
    double error9973()const{return 3 * stddev();}
    explicit NormalSummary(double theMean = 0, double theVariance = 0):
        mean(theMean), variance(theVariance){assert(variance >= 0);}
    NormalSummary operator-(NormalSummary const& b)const
        {return NormalSummary(mean - b.mean, variance + b.variance);}
    NormalSummary operator+=(NormalSummary const& b)
    {
        mean += b.mean;
        variance += b.variance;
        return *this;
    }
    NormalSummary operator*=(double a)
    {
        mean *= a;
        variance *= a * a;//change f code and description
        return *this;
    }
};

NormalSummary binomialRate(double k, int n)
{
    if(k < 0 || n < 1) return NormalSummary();
    double p = 1.0 * k/n;
    return NormalSummary(p, p * (1 - p)/n);
}
NormalSummary binomialCount(double k, int n)
    {return NormalSummary(k, k < 0 || n < 1 ? 0 : k * (1 - 1.0 * k/n));}

struct IncrementalStatistics
{
    double sum, squaredSum, minimum, maximum;
    long long n;
    IncrementalStatistics(): n(0), sum(0), squaredSum(0),
        minimum(numeric_limits<double>::max()), maximum(-minimum){}
    double getMean()const{return sum/n;}
    double getVariance()const{return n < 2 ? 0 :
        max(0.0, (squaredSum - sum * getMean())/(n - 1.0));}
    double stdev()const{return sqrt(getVariance());}
    void addValue(double x)
    {
        ++n;
        maximum = max(maximum, x);
        minimum = min(minimum, x);
        sum += x;
        squaredSum += x * x;
    }
    NormalSummary getSummary()
        {return NormalSummary(getMean(), getVariance()/n);}
    double error9973(){return getSummary().error9973();}
    double finiteSamplePlusMinusError01(double confidence = 0.9973)
    {//combined Hoeffding and empirical Bernstein
        assert(n > 1 && confidence > 0 && confidence < 1);
        double p = 1 - confidence, tempB = log(3/p);
        return min(sqrt(log(2/p)/2/n), (sqrt(2 * (n - 1) * getVariance() *
            tempB) + 3 * tempB)/n);
    }
    static double bernHelper(double x, double y, double sign)
    {
        return (x + 2 * y + sign * sqrt(x * (x + 4 * y * (1 - y))))/2/
            (1 + x);
    }
    pair<double, double> bernBounds(double confidence = 0.9973,
        int nFactor = 1)
    {
        assert(n > 0 && confidence > 0 && confidence < 1);
        double p = 1 - confidence, t = log(2/p)/n, yMin = getMean() - t/3,
            yMax = getMean() + t/3;
        return make_pair(yMin > 1 ? 0 : bernHelper(2 * t, yMin, -1),
            yMax > 1 ? 1 : bernHelper(2 * t, yMax, 1));
    }
    pair<double, double> combinedBounds(double confidence = 0.9973)
    {
        pair<double, double> result = bernBounds(confidence);
        double symBound = finiteSamplePlusMinusError01(confidence);
        result.first = max(result.first, getMean() - symBound);
        result.second = min(result.second, getMean() + symBound);
        return result;
    }
};

template<typename FUNCTION, typename STATISTIC> void MonteCarloSimulate(
    FUNCTION& f, STATISTIC& s, long long maxSimulations = 1000000,
    long long minSimulations = 1000,
    double precision = numeric_limits<double>::epsilon())
{
    for(int i = 0; i < maxSimulations && (i < minSimulations ||
        s.error9973() > precision); ++i) s.addValue(f());
}

template<typename FUNCTION,typename DATA> pair<double, pair<double, double> >
    bootstrap(Vector<DATA>const& data, FUNCTION const& f, int b = 10000,
    double confidence = 0.9973)
{
    assert(b > 2 && data.getSize() > 0);
    int tailSize = b * (1 - confidence)/2;
    if(tailSize < 1) tailSize = 1;
    if(tailSize > b/2 - 1) tailSize = b/2 - 1;
    Heap<double, ReverseComparator<double> > left;
    Heap<double> right;
    double sum = 0;
    for(int i = 0; i < b; ++i)
    {
        Vector<DATA> resample;
        for(int j = 0; j < data.getSize(); ++j)
            resample.append(data[GlobalRNG.mod(data.getSize())]);
        double value = f(resample);
        sum += value;
        if(left.getSize() < tailSize)
        {
            left.insert(value);
            right.insert(value);
        }
        else
        {
            if(value < left.getMin()) left.changeKey(0, value);
            else if(value > right.getMin()) right.changeKey(0, value);
        }
    }
    double Q = sum/b;
    return make_pair(Q, make_pair(Q - left.getMin(), right.getMin() - Q));
}

bool multipleComparison(Vector<NormalSummary>const& data,
    double meanPrecision = 0, double confidence = 0.9973)
{//smallest is best with precision meanPrecision
    assert(data.getSize() > 1);
    double probability = 1;
    for(int i = 1; i < data.getSize(); ++i)
    {
        NormalSummary diff = data[i] - data[0];
        probability -= (1 - approxNormalCDF(
            (diff.mean + meanPrecision)/sqrt(diff.variance)));
    }
    return probability >= confidence;
}

template<typename MULTI_FUNCTION> pair<Vector<IncrementalStatistics>, int>
    simulateSelectBest(MULTI_FUNCTION& f, int n0, int T,
    double meanPrecision = 0, double confidence = 0.9973)
{
    int D = f.getSize(), winner = -1;
    assert(D > 1 && n0 > 1 && T > n0 * D);
    Vector<IncrementalStatistics> data(D);
    for(int i = 0; i < D; ++i)
        for(int j = 0; j < n0; ++j) data[i].addValue(f(i));
    int k = n0 * D;
    for(; k < T; k += D)
    {
        Vector<NormalSummary> s;
        for(int i = 0; i < D; ++i) s.append(data[i].getSummary());
        int bestIndex = 0;
        double bestMean = s[0].mean;
        for(int i = 1; i < D; ++i)
            if(s[i].mean < bestMean) bestMean = s[bestIndex = i].mean;
        swap(s[0], s[bestIndex]);
        if(multipleComparison(s, meanPrecision, confidence))
        {
            winner = bestIndex;
            break;
        }
        for(int i = 0; i < D; ++i) data[i].addValue(f(i));
    }
    return make_pair(data, winner);
}

template<typename MULTI_FUNCTION> pair<Vector<IncrementalStatistics>, int>
    OCBA(MULTI_FUNCTION& f, int initialSims, int maxSims = 100000,
    double meanPrecision = 0, double confidence = 0.9973)
{
    int k = f.getSize(), winner = -1;
    assert(k > 1 && initialSims > 1 && maxSims > initialSims * k);
    Vector<IncrementalStatistics> data(k);
    for(int i = 0; i < k; ++i)
        for(int j = 0; j < initialSims; ++j) data[i].addValue(f(i));
    for(int j = initialSims * k; j < maxSims; ++j)
    {//put current best alternative in index 0
        Vector<NormalSummary> s;
        for(int i = 0; i < k; ++i) s.append(data[i].getSummary());
        int bestI = 0, bestRatioI = -1;
        double bestMean = s[0].mean, ratioSum = 0, bestRatio;
        for(int i = 1; i < k; ++i)
            if(s[i].mean < bestMean) bestMean = s[bestI = i].mean;
        swap(s[0], s[bestI]);
        //check if Bonferroni confidence is satisfied
        if(multipleComparison(s, meanPrecision, confidence))
        {
            winner = bestI;
            break;
        }
        //compute largest OCBA ratio
        for(int i = 1; i < k; ++i)
        {
            double meanDiff = s[i].mean - bestMean - meanPrecision, ratio =
                s[i].variance/(meanDiff * meanDiff);
            ratioSum += ratio * ratio/s[i].variance;
            if(bestRatioI == -1 || ratio > bestRatio)
            {
                bestRatio = ratio;
                bestRatioI = i;
            }
        }
        double ratioBest = sqrt(ratioSum * s[0].variance);
        if(ratioBest > bestRatio) bestRatioI = bestI;
        else if(bestRatioI == bestI) bestRatioI = 0;
        //simulate alterative with largest ratio
        data[bestRatioI].addValue(f(bestRatioI));
    }
    return make_pair(data, winner);
}

template<typename FUNCTION> struct SpeedTester
{
    FUNCTION f;
    SpeedTester(FUNCTION const& theFunction = FUNCTION()): f(theFunction){}
    int operator()()const
    {
        int now = clock();
        f();
        return clock() - now;
    }
};

unsigned char const SobolPolys[] = {0,1,1,2,1,4,2,4,7,11,13,14,1,13,16,19,22,
25,1,4,7,8,14,19,21,28,31,32,37,41,42,50,55,56,59,62,14,21,22,38,47,49,50,52,
56,67,70,84,97,103,115,122};
unsigned char const SobolDegs[] = {1,2,3,3,4,4,5,5,5,5,5,5,6,6,6,6,6,6,7,7,7,
7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8};
class Sobol
{//SobolPolys do not represent highest and lowest 1s
    enum{B = numeric_limits<double>::digits};
    unsigned long long k;
    Vector<unsigned long long> x, v;
    double factor;
    int index(int d, int b){return d * B + b;}
public:
    bool reachedLimit(){return k >= twoPower(B);}
    static int maxD(){return sizeof(SobolDegs);}
    Sobol(int d): factor(1.0/twoPower(B)), v(d * B, 0), x(d, 0), k(1)
    {
        assert(d <= maxD());
        for(int i = 0; i < d; ++i)
            for(int j = 0; j < B; ++j)
            {
                unsigned long long value;
                int l = j - SobolDegs[i];
                if(l < 0) value = (2 * j + 1) * twoPower(B - j - 1);
                else
                {
                    value = v[index(i, l)];
                    value ^= value/twoPower(SobolDegs[i]);
                    for(int k = 1; k < SobolDegs[i]; ++k)
                        if(Bits::get(SobolPolys[i], k - 1))
                            value ^= v[index(i, l + k)];
                }
                v[index(i, j)] = value;
            }
        next();
    }
    void next()
    {
        assert(!reachedLimit());
        for(int i = 0, c = rightmost0Count(k++); i < x.getSize(); ++i)
            x[i] ^= v[index(i, c)];
    }
    double getValue(int i){return x[i] * factor;}
};

template<typename PDF> class GridRWM
{
    PDF f;
    double x, fx, aFrom, aTo;
    int from, to;
    double sampleHelper(double a)
    {
        double xNew = x + GlobalRNG.uniform(-a, a), fxNew = f(xNew);
        if(fx * GlobalRNG.uniform01() <= fxNew)
        {
            x = xNew;
            fx = fxNew;
        }
        return x;
    }
public:
    GridRWM(double x0 = 0, PDF const& theF = PDF(), int from = -10,
        int to = 20): x(x0), f(theF), fx(f(x)), aFrom(pow(2, from)),
        aTo(pow(2, to)) {}
    double sample()
    {
        for(double a = aFrom; a < aTo; a *= 2) sampleHelper(a);
        return x;
    }
};

template<typename PDF> class MultidimGridRWM
{
    PDF f;
    Vector<double> x;
    double fx, aFrom, aTo, factor;
    Vector<double> sampleHelper(double a)
    {
        Vector<double> xNew = x;
        for(int i = 0; i < xNew.getSize(); ++i)
            xNew[i] += GlobalRNG.uniform(-a, a);
        double fxNew = f(xNew);
        if(fx * GlobalRNG.uniform01() <= fxNew)
        {
            x = xNew;
            fx = fxNew;
        }
        return x;
    }
public:
    MultidimGridRWM(Vector<double> const& x0, PDF const& theF = PDF(),
        int from = -10, int to = 20): x(x0), f(theF), fx(f(x)), aFrom(pow(2,
        from)), aTo(pow(2, to)), factor(pow(2, 1.0/x.getSize())) {}
    Vector<double> sample()
    {
        for(double a = aFrom; a < aTo; a *= factor) sampleHelper(a);
        return x;
    }
};

bool normalTestAreEqual(double m1, double v1, double m2, double v2,
    double z = 3)
{
    NormalSummary diff = NormalSummary(m1, v1) - NormalSummary(m2, v2);
    return abs(diff.mean) <= z * diff.stddev();
}

bool signTestAreEqual(double winCount1, double winCount2, double z = 3)
{
    double nGames = winCount1 + winCount2;
    return abs(winCount1 - nGames/2)/sqrt(nGames/4) < z;
}

struct SignedRankComparator
{
    typedef pair<double, double> P;
    int sign(P const& p)const{return p.first - p.second > 0 ? 1 : -1;}
    double diff(P const& p)const{return abs(p.first - p.second);}
    bool isLess(P const& lhs, P const& rhs)const
        {return diff(lhs) < diff(rhs);}
    bool isEqual(P const& lhs, P const& rhs)const
        {return diff(lhs) == diff(rhs);}
};
bool signedRankAreEqual(Vector<pair<double, double> > a, double z = 3)
{
    SignedRankComparator c;
    quickSort(a.getArray(), 0, a.getSize() - 1, c);
    int nP = a.getSize(), i = 0;
    //if odd number of 0's, drop first, distribute rest evenly
    while(i < a.getSize() && c.diff(a[i]) == 0) ++i;
    if(i % 2) --nP;
    double signedRankSum = 0, mean = nP * (nP + 1)/4;
    for(i = i % 2; i < a.getSize(); ++i)
    {//rank lookahead to scan for ties, then sum computation
        int j = i;
        while(i + 1 < a.getSize() && c.isEqual(a[i], a[i + 1])) ++i;
        double rank = (i + j)/2.0 + 1 + nP - a.getSize();
        while(j <= i) signedRankSum += c.sign(a[j++]) * rank;
    }
    return nP == 0 || (abs(signedRankSum) - mean)/
        sqrt(mean * (2 * nP + 1)/6) < z;
}

double evaluateChiSquaredCdf(double chi, int n)
{
    assert(chi >= 0 && n > 0);
    double m = 5.0/6 - 1.0/9/n - 7.0/648/n/n + 25.0/2187/n/n/n,
        q2 = 1.0/18/n + 1.0/162/n/n - 37.0/11664/n/n/n, temp = chi/n,
        x = pow(temp, 1.0/6) - pow(temp, 1.0/3)/2 + pow(temp, 1.0/2)/3;
    return approxNormalCDF((x - m)/sqrt(q2));
}
double chiSquaredP(Vector<int> const& counts,
    Vector<double> const& means, int degreesOfFreedom)
{
    double chiStat = 0;
    for(int i = 0; i < counts.getSize(); ++i)
        if(means[i] > 0 && counts[i] > 0) chiStat +=
            (counts[i] - means[i]) * (counts[i] - means[i])/means[i];
    return evaluateChiSquaredCdf(chiStat, degreesOfFreedom);
}

Vector<double> convertToRanks(Vector<double> a)
{//create index array, sort it, and convert indices into ranks
    int n = a.getSize();
    Vector<int> indices(n);
    for(int i = 0; i < n; ++i) indices[i] = i;
    IndexComparator<double> c(a.getArray());
    quickSort(indices.getArray(), 0, n - 1, c);
    for(int i = 0; i < n; ++i)
    {//rank lookahead to scan for ties, then change a entries
        int j = i;
        while(i + 1 < n && c.isEqual(i, i + 1)) ++i;
        double rank = (i + j)/2.0 + 1;
        for(; j <= i; ++j) a[indices[j]] = rank;
    }
    return a;
}

bool FriedmanAreAllEqual(Vector<Vector<double> > const& a,
    double conf = 0.95)
{//a[i] is vector of responses on domain i
    assert(a.getSize() > 0 && a[0].getSize() > 1);
    int n = a.getSize(), k = a[0].getSize();
    double totalRank = 0, sst = 0, sse = 0;
    Vector<double> optionTotalRanks(k);
    for(int i = 0; i < n; ++i)//calculate total and alternative rank sums
    {
        Vector<double> ri = convertToRanks(a[i]);
        for(int j = 0; j < k; ++j)
        {
            totalRank += ri[j];
            optionTotalRanks[j] += ri[j];
        }
    }
    totalRank /= n * k;
    optionTotalRanks *= 1.0/k/n;
    for(int i = 0; i < n; ++i)//calculate sums of squared ranks
    {
        sst += (optionTotalRanks[i] - totalRank) *
            (optionTotalRanks[i] - totalRank);
        Vector<double> ri = convertToRanks(a[i]);
        for(int j = 0; j < k; ++j)
            sse += (ri[j] - totalRank) * (ri[j] - totalRank);
    }
    sst *= n;
    sst /= n * (k - 1);
    return evaluateChiSquaredCdf(sst/sse, k - 1) > conf;
}

double evalFCdf(double x, int v1, int v2)
{//Paulson approximation
    assert(x >= 0 && v1 > 0 && v2 > 0);
    double temp1 = 2.0/v1/9, temp2 = 2.0/v2/9;
    return approxNormalCDF(((1 - temp2) * pow(x, 1.0/3) - (1 - temp1))/
        sqrt(temp2 * pow(x, 2.0/3) + temp1));
}
/*
double evalFCdf3(double x, int v1, int v2)
{//Scheffe-Tukey approximation
    double temp1 = 2 * v2, temp2 = v1 * x;
    double l = (temp1 + v1 - 2)/(temp1 + temp2);
    return evaluateChiSquaredCdf(temp2 * l, v1);
}

double evalFCdf2(double x, int v1, int v2)
{//Paulson approximation
    double temp1 = 2 * v2, temp2 = v1 * x/3;
    double l = (temp1 + temp2 + v1 - 2)/(temp1 + 4 * temp2);
    return evaluateChiSquaredCdf(v1 * x * l, v1);
}*/

double PearsonCorrelation(Vector<pair<double, double> > const& a)
{
    IncrementalStatistics x, y;
    for(int i = 0; i < a.getSize(); ++i)
    {
        x.addValue(a[i].first);
        y.addValue(a[i].second);
    }
    double covSum = 0;
    for(int i = 0; i < a.getSize(); ++i)
        covSum += (a[i].first - x.getMean()) * (a[i].second - y.getMean());
    return covSum/sqrt(x.getVariance() * y.getVariance());
}

double SpearmanCorrelation(Vector<pair<double, double> > a)
{
    Vector<double> x, y;
    for(int i = 0; i < a.getSize(); ++i)
    {
        x.append(a[i].first);
        y.append(a[i].second);
    }
    x = convertToRanks(x), y = convertToRanks(x);
    for(int i = 0; i < a.getSize(); ++i)
    {
        a[i].first = x[i];
        a[i].second = y[i];
    }
    return PearsonCorrelation(a);
}

}
#endif
