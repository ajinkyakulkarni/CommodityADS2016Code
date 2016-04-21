#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>
int evalCount = 0;
int evalCountG = 0;
#include "NumericalMethods.h"
#include "../Utils/DEBUG.h"
using namespace std;
using namespace igmdk;

struct SqrtFunction
{
    double value;
    SqrtFunction(double theValue):value(theValue){}
    double operator()(double x)const{return x * x - value;}
};

struct SqrtDerivative
{
    double operator()(double x)const{return 2 * x;}
};

double findSqrt(double x)
{
    return solveFor0(SqrtFunction(x), 1, x);
}

struct Stupid
{
    double x;
    double operator()()const
    {
        double d = (x - 100 + GlobalRNG.uniform01() * 10);
        return d * d;
    }
    Stupid(double theX):x(theX){}
};

struct PertQuad
{
    int D;
    typedef Vector<double> POINT;
    PertQuad(int theD):D(theD){}
    double operator()(POINT const& p)const
    {
        ++evalCount;
        double sum = 0;
        for(int i = 0; i < D; ++i)
        {
            sum += (i+1)*p[i]*p[i];
            double temp = 0;
            for(int j = 0; j < D; ++j) temp += p[j];
            sum += temp * temp/100;
        }
        return sum;
    }
    POINT start()const
    {
        POINT result(D, 0);
        for(int i = 0; i < D; ++i) result[i] = 0.5;
        return result;
    }
};

struct FH2
{
    int D;
    typedef Vector<double> POINT;
    FH2(int theD):D(theD){}
    double operator()(POINT const& p)const
    {
        ++evalCount;
        double sum = (p[0] - 3) * (p[0] - 3);
        for(int i = 1; i < D; ++i)
        {
            double temp = -1;
            for(int j = 0; j < i; ++j) temp += p[j];
            sum += temp * temp;
        }
        return sum;
    }
    POINT start()const
    {
        POINT result(D, 0);
        for(int i = 0; i < D; ++i) result[i] = 0.01;
        return result;
    }
};

struct Diagonal8
{
    int D;
    typedef Vector<double> POINT;
    Diagonal8(int theD):D(theD){}
    double operator()(POINT const& p)const
    {
        ++evalCount;
        double sum = 0;
        for(int i = 0; i < D; ++i)
        {
            sum += p[i]*(exp(p[i]) - 2 + p[i]);
        }
        return sum;
    }
    POINT start() const
    {
        POINT result(D, 0);
        for(int i = 1; i < D; ++i) result[i] = 1;
        return result;
    }
};

struct ExtTrid2
{
    int D;
    typedef Vector<double> POINT;
    ExtTrid2(int theD):D(theD){}
    double operator()(POINT const& p)const
    {
        ++evalCount;
        double sum = 0;
        for(int i = 0; i < D-1; ++i)
        {
            double temp = p[i] * p[i+1] - 1;
            sum += temp * temp + 0.1 * (p[i] + 1)*(p[i+1] + 1);
        }
        return sum;
    }
    POINT start()const
    {
        POINT result(D, 0);
        for(int i = 0; i < D; ++i) result[i] = 1;
        return result;
    }
};

struct ExtPenalty
{
    int D;
    typedef Vector<double> POINT;
    ExtPenalty(int theD):D(theD){}
    double operator()(POINT const& p)const
    {
        ++evalCount;
        double sum = 0, sum2 = 0;
        for(int i = 0; i < D; ++i)
        {
            sum2 += p[i] * p[i] - 0.25;
            if(i < D-1) sum += p[i] - 1;
        }
        return sum + sum2 * sum2;
    }
    POINT start()const
    {
        POINT result(D, 0);
        for(int i = 0; i < D; ++i) result[i] = i+1;
        return result;
    }
};

struct Fletcher
{
    int D;
    typedef Vector<double> POINT;
    Fletcher(int theD):D(theD){}
    double operator()(POINT const& p)const
    {
        ++evalCount;
        double sum = 0;
        for(int i = 1; i < D-1; ++i)
        {
            double temp = (p[i+1] - p[i] + 1 - p[i] * p[i]);
            sum += 100 * temp * temp;
        }
        return sum;
    }
    POINT start()const{return POINT(D, 0);}
};

struct Staircase2
{
    int D;
    typedef Vector<double> POINT;
    Staircase2(int theD):D(theD){}
    POINT start()const{return POINT(D, 0);}
    double operator()(POINT const& p)const
    {
        ++evalCount;
        double sum = 0;
        for(int i = 0; i < D; ++i)
        {
            double sum2 = -i;
            for(int j = 0; j <= i; ++j)
            {
                sum2 += p[j];
            }
            sum += sum2 * sum2;
        }
        return sum;
    }
};

struct BDQRTIC
{
    int D;
    typedef Vector<double> POINT;
    BDQRTIC(int theD):D(theD){}
    double operator()(POINT const& p)const
    {
        ++evalCount;
        assert(D > 4);
        double sum = 0;
        for(int i = 0; i < D-4; ++i)
        {
            double temp = -4 * p[i] + 3,
                temp2 = 5 * p[D-1] * p[D-1];
            for(int j = 0; j < 4; ++j) temp2 += (i+1) * p[i+j] * p[i+j];
            sum += temp * temp + temp2 * temp2;
        }
        return sum;
    }
    POINT start()const
    {
        POINT result(D, 0);
        for(int i = 0; i < D; ++i) result[i] = 1;
        return result;
    }
};

struct DQDRTIC
{
    int D;
    typedef Vector<double> POINT;
    DQDRTIC(int theD):D(theD){}
    double operator()(POINT const& p)const
    {
        ++evalCount;
        assert(D > 2);
        double sum = 0;
        for(int i = 0; i < D-2; ++i)
        {
            sum += p[i] * p[i] + 100 * p[i+1] * p[i+1] + 100 * p[i+2] * p[i+2];
        }
        return sum;
    }
    POINT start()const
    {
        POINT result(D, 0);
        for(int i = 0; i < D; ++i) result[i] = 3;
        return result;
    }
};

struct EDENSCH
{
    int D;
    typedef Vector<double> POINT;
    EDENSCH(int theD):D(theD){}
    double operator()(POINT const& p)const
    {
        ++evalCount;
        assert(D > 2);
        double sum = 16;
        for(int i = 0; i < D-1; ++i)
        {
            double temp = p[i] - 2, temp2 = p[i] * p[i+1] - 2 * p[i+1], temp3 = p[i+1] + 1;
            sum += temp * temp + temp2 * temp2 + temp3 * temp3;
        }
        return sum;
    }
    POINT start()const
    {
        POINT result(D, 0);
        for(int i = 0; i < D; ++i) result[i] = 0;
        return result;
    }
};

struct Stupid3
{
    typedef Point<double, 2> POINT;
    static POINT start(){return POINT(400, 200);}
    double operator()(Point<double, 2>const& p)const
    {
        ++evalCount;
        return sqrt(p[0] * p[0] + p[1] * p[1]) +
            sqrt((p[0] - 300) * (p[0] - 300) + (p[1] - 400) * (p[1] - 400)) +
            sqrt((p[0] - 700) * (p[0] - 700) + (p[1] - 300) * (p[1] - 300));
    }

};

struct Stupid3Grad
{
    Point<double, 2> operator()(Point<double, 2>const& p)const
    {
        ++evalCountG;
        assert(p[0] != 0 && p[0] != 300 && p[0] != 700 && p[1] != 0 && p[1] != 300 && p[1] != 400);
        double temp1 = sqrt(p[0] * p[0] + p[1] * p[1]), temp2 =
            sqrt((p[0] - 300) * (p[0] - 300) + (p[1] - 400) * (p[1] - 400)),
            temp3 = sqrt((p[0] - 700) * (p[0] - 700) + (p[1] - 300) * (p[1] - 300));
        return Point<double, 2>(p[0]/temp1 + (p[0] - 300)/temp2 + (p[0] - 700)/temp3,
            p[1]/temp1 + (p[1] - 400)/temp2 + (p[1] - 300)/temp3);
    }
};

struct Killer
{
    typedef Point<double, 2> POINT;
    static POINT start(){return POINT(-1, -1);}
    double operator()(Point<double, 2>const& p)const
    {
        ++evalCount;
        return max((p[0]-1)*(p[0]-1)+(p[1]+1)*(p[1]+1),(p[0]+1)*(p[0]+1)+(p[1]-1)*(p[1]-1));
    }
};

struct StupidD
{
    int D;
    typedef Vector<double> POINT;
    POINT start()const{return POINT(D, 0);}
    StupidD(int theD):D(theD){}
    double operator()(POINT const& p)const
    {
        ++evalCount;
        double sum = 200;
        for(int i = 0; i < D; ++i)
            sum += (p[i] - 400.0/3) * (p[i] - 400.0/3);
        return  sum;
    }
};

struct Log30
{
    int D;
    typedef Vector<double> POINT;
    POINT start()const{return POINT(D, 0);}
    Log30(int theD):D(theD){}
    double operator()(POINT const& p)const
    {
        ++evalCount;
        double sum = 0;
        for(int i = 0; i < D; ++i)
        {
            double temp = max(1.0, log(p[i]+10)) - 30;
            sum += temp * temp;
        }
        return  sum;
    }
};

template<int D> struct StupidY
{
    typedef Point<double, D> POINT;
    POINT p;
    double operator()()const
    {
        double sum = 200;
        for(int i = 0; i < D; ++i)
        {
            double temp = p[i] - 400.0/3 + GlobalRNG.uniform01();
            sum += temp * temp;
        }
        return  sum;
    }
    StupidY(POINT const& theP):p(theP)
    {

    }
};

template<int D> struct StupidZ
{
    typedef Point<double, D> POINT;
    static POINT start(){return POINT();}
    double operator()(Point<double, D>const& p)const
    {
        ++evalCount;
        StupidY<D> f(p);
        IncrementalStatistics s;
        MonteCarloSimulate(f, s, 1000, 1000);
        return s.getSummary().mean();
    }
};

struct Stupid2Grad
{
    Point<double, 2> operator()(Point<double, 2>const& p)const
    {
        ++evalCountG;
        return Point<double, 2>(2*(p[0] - 400.0/3),
            2 * (p[1] - 400.0/3));
    }
};

struct RosenbrockGrad
{
    Point<double, 2> operator()(Point<double, 2>const& p)const
    {
        ++evalCountG;
        double temp = p[1] - p[0] * p[0];
        return Point<double, 2>(-2*(1 - p[0]) - 400.0 * p[0] * temp, 200 * temp);
    }
};

struct Rosenbrock
{
    int D;
    typedef Vector<double> POINT;
    Rosenbrock(int theD):D(theD){}
    POINT start()const{return POINT(D, 0);}
    double operator()(POINT const& p)const
    {
        ++evalCount;

        double sum = 0;
        for(int i = 1; i < D; ++i)
        {
            double temp1 = 1 - p[i-1], temp2 = p[i] - p[i-1] * p[i-1];
            sum += temp1*temp1 + 100 * temp2 * temp2;
        }
        return sum;
    }
};

template<typename POINT> void debugResult(pair<POINT, double> const& result)
{
    //for(int i = 0; i < result.first.getSize(); ++i) DEBUG(result.first[i]);
    DEBUG(result.second);
    DEBUG(evalCount);
    evalCount = 0;
    DEBUG(evalCountG);
    evalCountG = 0;
}

template<typename TESTCASE> void testAllSolvers(TESTCASE const& f)
{
    typedef typename TESTCASE::POINT POINT;
    typedef GradientFunctor<TESTCASE, POINT> GRADIENT;
    GRADIENT g(f);
    /*DEBUG("metaSPSA");
    debugResult(metaSPSA(f.start(), f));*/
    /*GRADIENT g(f);
    DEBUG("SPSA");
    debugResult(SPSA(f.start(), f, 1000000, 0.001));*/
    /*DEBUG("DSRS");
    debugResult(DSRS(f.start(), f));*/

    /*DEBUG("ILSNelderMead");
    ILSNelderMead<POINT, TESTCASE> ilsnm;
    debugResult(ilsnm.minimize(TESTCASE::start()));*/
    /*DEBUG("NelderMead");
    NelderMead<POINT, TESTCASE> nm(f.start().getSize(), f);
    debugResult(nm.minimize(f.start(), 100000, numeric_limits<double>::epsilon(), f.start().getSize()));*/
    /*DEBUG("RestartedNelderMead");
    NelderMead<POINT, TESTCASE> nm2(f.start().getSize(), f);
    debugResult(nm2.restartedMinimize(f.start(), 100000, numeric_limits<double>::epsilon(), 100, f.start().getSize()));*/
    DEBUG("LBFGSMinimize");
    debugResult(LBFGSMinimize(f.start(), f, g));
}


struct FunctionOne
{
    template<typename POINT> double operator()(POINT const& point)const{return 1;}
};

void testIntegrate()
{
    double result = Trapezoid<SqrtFunction>::integrate(SqrtFunction(0), 0, 1);
    DEBUG(result);
}

struct UnitSquareTest
{
    template<typename POINT> double operator()(POINT const& point)const
    {return point[0] * point[0] + point[1] * point[1] <= 1;}
};

void testMCI()
{
    Point<double, 2> a(-1, -1), b(1, 1);
    int n = 10000000;
    //sobol 3.1415 va 2.7-e6
    //random 3.14463 va 2.7-e6
    double value = MonteCarloIntegrate<Point<double, 2>, UnitSquareTest, FunctionOne>(make_pair(a, b), n);
    DEBUG(value);
}

struct ABS
{
    double operator()(double x)const
    {
        return abs(x);
    }
};

struct ABS2
{
    double operator()(double x)const
    {
        return abs(x-1);
    }
};

void testAllFunctions(int D)
{
    DEBUG(D);
    DEBUG("Stupid");
    testAllSolvers(StupidD(D));
    DEBUG("Fletcher");
    testAllSolvers(Fletcher(D));
    DEBUG("Rosenbrock");
    testAllSolvers(Rosenbrock(D));
    DEBUG("ExtTrid2");
    testAllSolvers(ExtTrid2(D));
    DEBUG("Diagonal8");
    testAllSolvers(Diagonal8(D));
    DEBUG("ExtPenalty");
    testAllSolvers(ExtPenalty(D));
    DEBUG("BDQRTIC");
    //DEBUG("N/A");//for D = 2
    testAllSolvers(BDQRTIC(D));
    DEBUG("DQDRTIC");
    //DEBUG("N/A");//for D = 2
    testAllSolvers(DQDRTIC(D));
    DEBUG("EDENSCH");
    //DEBUG("N/A");//for D = 2
    testAllSolvers(EDENSCH(D));
    DEBUG("Staircase2");
    testAllSolvers(Staircase2(D));
    DEBUG("PertQuad");
    testAllSolvers(PertQuad(D));
    DEBUG("FH2");
    testAllSolvers(FH2(D));
    /*DEBUG("Log30");
    testAllSolvers(Log30(D));*/
}

struct EPX
{
    double operator()(double x)const{return exp(x);}
};

void testChebApprox()
{
    EPX z;
    ChebFunction cf(z, 16);
    DEBUG(cf(0.3));
    DEBUG(z(0.3));
    ChebFunction ci = cf.integral();

    DEBUG(ci(-1));
    DEBUG(ci(1));
    DEBUG(ci(1) - ci(-1));
    DEBUG(cf.integrateM11());
    Poly<double> p = cf.poly();
    p.debug();
    DEBUG(p(0.3));
    DEBUG(Trapezoid<EPX>::integrate(z, -1, 1));
    Poly<double> ip = p.integrate();
    DEBUG(ip(1) - ip(-1));
    DEBUG(z(1) - z(-1));
    ChebFunction di = ci.derivative();
    DEBUG(di(0.3));
    DEBUG(integrate(z, -1, 1));
}

int main()
{
    /*DEBUG(findSqrt(3));

    testMCI();
    testIntegrate();
    DEBUG(numeric_limits<double>::min());
    DEBUG(numeric_limits<double>::max());
    DEBUG(numeric_limits<double>::epsilon());*/

    /*DEBUG("Killer");
    testAllSolvers(Killer());
    DEBUG("StupidZ<2>");
    testAllSolvers(StupidZ<2>());
    testAllFunctions(2);*/
    //testAllFunctions(5);
    //testAllFunctions(10);
    //testAllFunctions(100);
    //testAllFunctions(1000);
    //testAllFunctions(10000);

    /*Killer k;
    DEBUG(k(Point<double, 2>(-1, -1)));
    DEBUG(k(Point<double, 2>(-1, 0)));
    DEBUG(k(Point<double, 2>(0, -1)));*/

    /*typedef Rosenbrock<2> TESTCASE;
    typedef TESTCASE::POINT POINT;
    typedef RosenbrockGrad GRADIENT;
    GRADIENT g;*/
    /*DEBUG("conjugateGradientMinimize");
    debugResult(conjugateGradientMinimize(TESTCASE::start(), TESTCASE(), g));*/
    testChebApprox();
    return 0;
}
