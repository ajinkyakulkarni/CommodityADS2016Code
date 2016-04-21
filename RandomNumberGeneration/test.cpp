#include "Random.h"
#include "Statistics.h"
#include "../Utils/Debug.h"
#include "../NumericalMethods/Matrix.h"
using namespace igmdk;

struct meder
{
    double operator()(Vector<double> observations)const
    {
        quickSort(observations.getArray(), 0, observations.getSize());
        return observations[observations.getSize()/2];
    }
};

struct meaner
{
    double operator()(Vector<double> const& observations)const
    {
        double result = observations[0];
        for(int i = 1; i < observations.getSize(); ++i)
        {
            result += observations[i];
        }
        return result / observations.getSize();
    }
};

struct maxer
{
    double operator()(Vector<double> const& observations)const
    {
        double result = observations[0];
        for(int i = 1; i < observations.getSize(); ++i)
        {
            if(observations[i] > result) result = observations[i];
        }
        return result;
    }
};

void testSumHeap()
{
    int N = 1000000;
    SumHeap<double> sh;
    sh.add(0);
	sh.add(1.0/36);
	sh.add(2.0/36);
	sh.add(3.0/36);
	sh.add(4.0/36);
	sh.add(5.0/36);
	sh.add(6.0/36);
	sh.add(5.0/36);
	sh.add(4.0/36);
	sh.add(3.0/36);
	sh.add(2.0/36);
	sh.add(1.0/36);
    int sum = 0;
    for(int i = 0 ; i < N; ++i) sum += sh.next();
    DEBUG(sum*1.0/N);
}

struct XYZ
{
    double operator()()const
    {
        return GlobalRNG.bernoulli(0.95);
        //return GlobalRNG.uniform01();
    }
};
struct OCBATest
{
    int getSize(){return 6;}
    double operator()(int i)
    {
        if(i == 0) return GlobalRNG.normal(0, 10);
        if(i == 1) return GlobalRNG.normal(0.1, 9);
        if(i == 2) return GlobalRNG.normal(0.2, 8);
        if(i == 3) return GlobalRNG.normal(0.3, 7);
        if(i == 4) return GlobalRNG.normal(0.4, 6);
        else return GlobalRNG.normal(0.5, 5);
    }
};

void testOCBA()
{
    OCBATest t;
    DEBUG(OCBA(t, 1000, 1000000).second);
    DEBUG(simulateSelectBest(t, 1000, 1000000).second);
}

struct Normal01SemiPDF
{
    double mean, variance;
    Normal01SemiPDF(double theMean = 100, double theVariance = 10000):
        mean(theMean), variance(theVariance){}
    double operator()(double x)const{x -= mean; return exp(-x * x/2/variance);}
};

struct MultivarNormalSemiPDF
{
    double operator()(Vector<double> x)const
    {
        Matrix<double> b2(3, 3);
        b2(0, 0) = 1;
        b2(0, 1) = 4;
        b2(0, 2) = 5;
        b2(1, 0) = 4;
        b2(1, 1) = 20;
        b2(1, 2) = 32;
        b2(2, 0) = 5;
        b2(2, 1) = 32;
        b2(2, 2) = 64;
        LUP<double> lup(b2);
        Matrix<double> inv = lup.inverse();
        //inv.debug();
        //for(int i = 0; i < x.getSize(); ++i) x[i] -= 100;
        //for(int i = 0; i < x.getSize(); ++i) DEBUG(x[i]);
        //DEBUG(inv * x * x/(-2));
        //DEBUG(exp(inv * x * x/(-2)));
        //system("PAUSE");
        return exp(inv * x * x/(-2));
    }
};

void testVectorMagicMCMC2()
{
    MultidimGridRWM<MultivarNormalSemiPDF> g(Vector<double>(3, 0.5));
    int n = 10000;
    for(int i = 0; i < n; ++i) g.sample();

    Vector<double> sum(3, 0);
    Matrix<double> outerSum(3, 3);
    for(int i = 0; i < n; ++i)
    {
        Vector<double> x = g.sample();
        sum += x;
        outerSum += outerProduct(x, x);
    }
    Vector<double> mean = sum * (1.0/n);
    for(int i = 0; i < 3; ++i) DEBUG(mean[i]);
    Matrix<double> cov = (outerSum - outerProduct(mean, sum)) * (1.0/(n - 1));
    cov.debug();
}

void testMCMCMagic()
{
    GridRWM<Normal01SemiPDF> s;
    int n = 10;
    for(int i = 0; i < n; ++i) s.sample();
    IncrementalStatistics z;

    for(int i = 0; i < n; ++i) z.addValue(s.sample());
    DEBUG(z.getMean());
    DEBUG(z.getVariance());

}

void testStatTests()
{//should be 0.95 for all
    DEBUG(evaluateChiSquaredCdf(3.84, 1));
    DEBUG(evaluateChiSquaredCdf(11.1, 5));
    DEBUG(evaluateChiSquaredCdf(18.3, 10));
    DEBUG(evaluateChiSquaredCdf(31.4, 20));
    DEBUG(evaluateChiSquaredCdf(124, 100));
}

void testF()
{//should be 0.95 for all
    DEBUG(evalFCdf(161.4, 1, 1));
    DEBUG(evalFCdf(19, 2, 2));
    DEBUG(evalFCdf(9.28, 3, 3));
    DEBUG(evalFCdf(6.39, 4, 4));
    DEBUG(evalFCdf(5.05, 5, 5));
    DEBUG(evalFCdf(2.98, 10, 10));
    DEBUG(evalFCdf(2.12, 20, 20));
}
/*
void testF2()
{//should be 0.95 for all
    DEBUG(evalFCdf2(161.4, 1, 1));
    DEBUG(evalFCdf2(19, 2, 2));
    DEBUG(evalFCdf2(9.28, 3, 3));
    DEBUG(evalFCdf2(6.39, 4, 4));
    DEBUG(evalFCdf2(5.05, 5, 5));
    DEBUG(evalFCdf2(2.98, 10, 10));
    DEBUG(evalFCdf2(2.12, 20, 20));
}*/

void testWilcoxon()
{
    Vector<pair<double, double> > a;
    a.append(pair<double, double>(125, 110));
    a.append(pair<double, double>(115, 122));
    a.append(pair<double, double>(130, 125));
    a.append(pair<double, double>(140, 120));
    a.append(pair<double, double>(140, 140));
    a.append(pair<double, double>(115, 124));
    a.append(pair<double, double>(140, 123));
    a.append(pair<double, double>(125, 137));
    a.append(pair<double, double>(140, 135));
    a.append(pair<double, double>(135, 145));
    signedRankAreEqual(a, 2);
}

void DDDAlias()
{
    Vector<double> probabilities;
    for(int i = 0; i < 5; ++i) probabilities.append((i-2)*(i-2)+1);
    normalizeProbs(probabilities);
    AliasMethod alias(probabilities);
    cout << "breakpoint" << endl;
}

void DDDSumHeap()
{
    Vector<double> probabilities;
    for(int i = 0; i < 5; ++i) probabilities.append((i-2)*(i-2)+1);
    normalizeProbs(probabilities);
    SumHeap<double> sumHeap;
    for(int i = 0; i < 5; ++i) sumHeap.add(probabilities[i]);
    cout << "breakpoint" << endl;
}

int main(int argc, char *argv[])
{
    DDDAlias();
    DDDSumHeap();
    testWilcoxon();
    return 0;
    testF();
    testStatTests();
    //return 0;
    DEBUG(signTestAreEqual(7, 10, 0.5));
    testSumHeap();
    //return 0;
    XYZ xyz;
    IncrementalStatistics si;
    MonteCarloSimulate(xyz, si, 1000);
    NormalSummary s = si.getSummary();
    DEBUG(s.mean);
    DEBUG(s.error9973());
    DEBUG(si.finiteSamplePlusMinusError01());
    pair<double, double> bb = si.bernBounds();
    DEBUG(si.getMean() - bb.first);
    DEBUG(bb.second - si.getMean());
    bb = si.combinedBounds();
    DEBUG(si.getMean() - bb.first);
    DEBUG(bb.second - si.getMean());
    return 0;
    Vector<double> data;
    //data.append(1000);
    //data.append(2000);
    meder med;
    for(int i = 0; i < 1000; ++i) data.append(GlobalRNG.uniform01());
    pair<double, pair<double, double> > r2 = bootstrap(data, med, 10000);
    double medd = med(data);
    DEBUG(medd);
    DEBUG(r2.first);
    DEBUG(r2.second.first);
    DEBUG(r2.second.second);


    DEBUG(Sobol::maxD());
    DEBUG(sizeof(SobolPolys));
    DEBUG(sizeof(SobolDegs));
    Sobol so(1);
    for(int i = 0; i < 10; ++i)
    {
        for(int j = 0; j < 1; ++j)
            DEBUG(so.getValue(j));
        so.next();
    }
    //ARC4 x;
    MRG32k3a x;
    x.jumpAhead();
    unsigned long long N = 1 << 3;
    unsigned long long dummy = 0;
    while(N--) dummy += x.next();
    DEBUG(dummy);

    SumHeap<double> st;
    st.add(0.2);
    st.add(0.2);
    st.add(0.2);
    st.add(0.2);
    st.add(0.2);
    DEBUG(st.total());
    DEBUG(st.find(0.1));
    DEBUG(st.find(0.3));
    DEBUG(st.find(0.6));
    DEBUG(st.find(0.9));

    DEBUG(st.find(0));
    DEBUG(st.find(0));
    DEBUG(st.find(0.5));
    DEBUG(st.find(1));

    DEBUG(st.cumulative(0));
    DEBUG(st.cumulative(1));
    DEBUG(st.cumulative(2));
    DEBUG(st.cumulative(3));
    DEBUG(st.cumulative(4));
    DEBUG(st.find(st.cumulative(0)));
    DEBUG(st.find(st.cumulative(1)));
    DEBUG(st.find(st.cumulative(2)));
    DEBUG(st.find(st.cumulative(3)));
    DEBUG(st.find(st.cumulative(4)));
    for(int i = 0; i < 100; ++i)
    {
        DEBUG(st.next());
    }

    clock_t start = clock();
    //Xorshift random;
    //Xorshift64 random;
    //ImprovedXorshift random;
    QualityXorshift64 random;
    //ImprovedXorshift64 random;
    unsigned long long sum = 0;
    for(int i = 0; i < 1000000; ++i)
    {
		sum += random.next();
	}
	DEBUG(sum);
	clock_t end = clock();
	int time = (end - start);
    cout << "IX: " << time << endl;

    if(true)
    {
        DEBUG(GlobalRNG.uniform01());
        DEBUG(GlobalRNG.uniform(10, 20));
        DEBUG(GlobalRNG.normal01());
        DEBUG(GlobalRNG.normal(10, 20));
        DEBUG(GlobalRNG.exponential(1));
        DEBUG(GlobalRNG.gamma1(0.5));
        DEBUG(GlobalRNG.gamma1(1.5));
        DEBUG(GlobalRNG.weibull1(20));
        DEBUG(GlobalRNG.erlang(10, 2));
        DEBUG(GlobalRNG.chiSquared(10));
        DEBUG(GlobalRNG.t(10));
        DEBUG(GlobalRNG.logNormal(10, 20));
        DEBUG(GlobalRNG.beta(0.5, 0.5));
        DEBUG(GlobalRNG.F(10 ,20));
        DEBUG(GlobalRNG.cauchy(0, 1));
        DEBUG(GlobalRNG.binomial(0.7, 20));
        DEBUG(GlobalRNG.geometric(0.7));
        DEBUG(GlobalRNG.poisson(0.7));
		//system("PAUSE");
	}
	int M = 100000;
	double average = 0;
	for(int i = 0; i < M; ++i)
	{
        average += GlobalRNG.beta(0.5, 0.5);
	}
	DEBUG(average/M);
	DEBUG(approxNormal2SidedConf(3));
	Vector<NormalSummary> hha;
	NormalSummary n1(10, 100.0/8);
	NormalSummary n2(20, 81.0/8);
	NormalSummary n3(22, 144.0/8);
	hha.append(n1);
	hha.append(n2);
	hha.append(n3);
	DEBUG(multipleComparison(hha, 0, 0.95));
	testOCBA();
	Vector<NormalSummary> hha2;
	NormalSummary n4(1, 0.1);
	NormalSummary n5(2, 0.2);
	NormalSummary n6(3, 0.3);
	hha2.append(n4);
	hha2.append(n5);
	hha2.append(n6);
	DEBUG(multipleComparison(hha2, 0, 0.95));
    testVectorMagicMCMC2();
    testMCMCMagic();
    unsigned long long zzzz = 12345678ull * 2685821657736338717ull;
	DEBUG(12345678ull * 2685821657736338717ull);
	DEBUG(zzzz * 5.42101086242752217E-20);
    return 0;
}
