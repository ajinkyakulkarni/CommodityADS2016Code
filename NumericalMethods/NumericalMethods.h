#ifndef NUMERICAL_METHODS_H
#define NUMERICAL_METHODS_H
#include <cmath>
#include <complex>
#include <iomanip>
#include "../Utils/Bits.h"
#include "../Utils/Vector.h"
#include "../Utils/Utils.h"
#include "../Utils/Sort.h"
#include "../Utils/Queue.h"
#include "../RandomNumberGeneration/Statistics.h"
#include "../RandomTreap/Treap.h"
#include "../ComputationalGeometry/Point.h"
#include "../Optimization/Metaheuristics.h"
#include "Matrix.h"
namespace igmdk{

double globalSafeDelta = sqrt(numeric_limits<double>::epsilon());

bool isELess(double a, double b,
    double eRelAbs = numeric_limits<double>::epsilon())
    {return a < b && b - a >= eRelAbs * (abs(a) + abs(b) + 1);}
bool isEEqual(double a, double b,
    double eRelAbs = numeric_limits<double>::epsilon())
    {return !isELess(a, b, eRelAbs) && !isELess(b, a, eRelAbs);}

bool haveDifferentSign(double a, double b){return (a < 0) != (b < 0);}
template<typename FUNCTION> double solveFor0(FUNCTION const& f, double xLeft,
    double xRight)
{
    double yLeft = f(xLeft), xMiddle;
    assert(xRight >= xLeft && haveDifferentSign(yLeft, f(xRight)));
    for(;;)
    {
        xMiddle = (xLeft + xRight)/2;
        double yMiddle = f(xMiddle), prevDiff = xRight - xLeft;
        if(haveDifferentSign(yLeft, yMiddle)) xRight = xMiddle;
        else
        {
            xLeft = xMiddle;
            yLeft = yMiddle;
        }
        if(xRight - xLeft >= prevDiff) break;
    }
    return xMiddle;
}

template<typename FUNCTION> double estimateDerivativeCD(double x,
    FUNCTION const& f, double relAbsDelta = globalSafeDelta)
{
    double h = globalSafeDelta * max(1.0, abs(x));
    return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h))/h/12;
}

template<typename FUNCTION, typename POINT> POINT estimateGradientFD(POINT x,
    FUNCTION const& f, double relAbsDelta = globalSafeDelta)
{
    POINT result = x;
    double y = f(x);
    for(int i = 0; i < x.getSize(); ++i)
    {
        double h = relAbsDelta * max(1.0, abs(x[i]));
        x[i] += h;
        result[i] = (f(x) - y)/h;
        x[i] -= h;
    }
    return result;
}

template<typename FUNCTION, typename POINT> struct GradientFunctor
{
    FUNCTION f;
    double delta;
    GradientFunctor(FUNCTION const& theF,
        double relAbsDelta = globalSafeDelta): f(theF), delta(relAbsDelta) {}
    POINT operator()(POINT const& p)const
        {return estimateGradientFD(p, f, delta);}
};

template<typename FUNCTION, typename DERIVATIVE> double solveFor0Newton(
    FUNCTION const& f, DERIVATIVE const& derivative, double x)
{
    for(double prevDiff = numeric_limits<double>::max();;)
    {
        double y = f(x), oldX = x;
        x -= y/derivative(x);
        double diff = abs(x - oldX);
        if(diff >= prevDiff) break;
        prevDiff = diff;
    }
    return x;
}

template<typename FUNCTION> pair<double, double> minimizeGS(
    FUNCTION const& f, double xLeft, double xRight,
    double relAbsXPrecision = globalSafeDelta)
{
    assert(isfinite(xLeft) && isfinite(xRight) && xLeft <= xRight &&
        relAbsXPrecision >= numeric_limits<double>::epsilon());
    double GR = 0.618, xMiddle = xLeft * GR + xRight * (1 - GR),
        yMiddle = f(xMiddle);
    while(isELess(xLeft, xRight, relAbsXPrecision))
    {
        bool chooseR = xRight - xMiddle > xMiddle - xLeft;
        double xNew = GR * xMiddle + (1 - GR) *
            (chooseR ? xRight : xLeft), yNew = f(xNew);
        if(yNew < yMiddle)
        {
            (chooseR ? xLeft : xRight) = xMiddle;
            xMiddle = xNew;
            yMiddle = yNew;
        }
        else (chooseR ? xRight : xLeft) = xNew;
    }
    return make_pair(xMiddle, yMiddle);
}
int roundToNearestInt(double x)
{
    return int(x > 0 ? x + 0.5 : x - 0.5);
}
template<typename FUNCTION> pair<int, double> minimizeGSDiscrete(
    FUNCTION const& f, int xLeft, int xRight)
{
    assert(isfinite(xLeft) && isfinite(xRight) && xLeft <= xRight);
    double GR = 0.618;
    int xMiddle = roundToNearestInt(xLeft * GR + xRight * (1 - GR));
    double yMiddle = f(xMiddle);
    while(xLeft < xRight)
    {
        bool chooseR = xRight - xMiddle > xMiddle - xLeft;
        int xNew = roundToNearestInt(GR * xMiddle + (1 - GR) *
            (chooseR ? xRight : xLeft));
        double yNew = xNew == xMiddle ? yMiddle : f(xNew);
        if(yNew < yMiddle)
        {
            (chooseR ? xLeft : xRight) = xMiddle;
            xMiddle = xNew;
            yMiddle = yNew;
        }
        else (chooseR ? xRight : xLeft) = xNew;
    }
    return make_pair(xMiddle, yMiddle);
}

template<typename FUNCTION> double findIntervalBound(FUNCTION const& f,
    double guess, double d)
{//run with d < 0 for the left bound and d > 0 for the right bound
    for(double yBest = f(guess); d * 2 != d; d *= 2)//handle d = 0 and inf
    {
        //DEBUG(guess);
        //DEBUG(yBest);
        //DEBUG(guess + d);
        double yNext = f(guess + d);
        //DEBUG(yNext);
        if(yNext >= yBest) break;
        yBest = yNext;
    }
    return guess + d;
}

template<typename POINT, typename FUNCTION> class NelderMead
{
    FUNCTION f;
    int D;
    POINT vertexSum;//incremental centroid
    typedef pair<POINT, double> P;
    Vector<P> simplex;
    double scale(P& high, double factor)
    {
        P result = high;
        //affine combination of the high point and the
        //centroid of the remaining vertices
        //centroid = (vertexSum - high)/D and
        //result = centroid * (1 - factor) + high * factor
        double centroidFactor = (1 - factor)/D;
        result.first = vertexSum * centroidFactor +
            high.first * (factor - centroidFactor);
        result.second = f(result.first);
        if(result.second < high.second)
        {//accept scaling if improving
            vertexSum += result.first - high.first;
            high = result;
        }
        return result.second;
    }
public:
    NelderMead(int theD, FUNCTION const& theFunction = FUNCTION()):
        D(theD), f(theFunction), simplex(D + 1, D + 1){}

    P minimize(POINT const& initialGuess, int maxIterations = 10000,
        double yPrecision = globalSafeDelta, double step = 1)
    {
        vertexSum = initialGuess;
        for(int i = 0; i < D; ++i) vertexSum[i] = 0;
        for(int i = 0; i <= D; ++i)
        {
            simplex[i].first = initialGuess;
            if(i > 0)simplex[i].first[i - 1] += GlobalRNG.uniform01() * step;
            simplex[i].second = f(simplex[i].first);
            vertexSum += simplex[i].first;
        }
        for(;;)
        {//calculate high, low, and nextHigh, which must be all different
            int high = 0, nextHigh = 1, low = 2;
            if(simplex[high].second < simplex[nextHigh].second)
                swap(high, nextHigh);
            if(simplex[nextHigh].second < simplex[low].second)
            {
                swap(low, nextHigh);
                if(simplex[high].second < simplex[nextHigh].second)
                    swap(high, nextHigh);
            }
            for(int i = 3; i <= D; ++i)
            {
                if(simplex[i].second < simplex[low].second) low = i;
                else if(simplex[i].second > simplex[high].second)
                {
                    nextHigh = high;
                    high = i;
                }
                else if(simplex[i].second > simplex[nextHigh].second)
                    nextHigh = i;
            }
            if(!maxIterations-- || !isELess(simplex[low].second,
                simplex[high].second, yPrecision)) return simplex[low];
            //try to reflect
            double value = scale(simplex[high], -1);
            //try to double if better than low
            if(value <= simplex[low].second) scale(simplex[high], 2);
            else if(value >= simplex[nextHigh].second)
            {//try reflected/unrefrected halving if accepted/rejected value
                double yHi = simplex[high].second;
                if(scale(simplex[high], 0.5) >= yHi)
                {//contract all to get rid of the high point
                    vertexSum = simplex[low].first;
                    for(int i = 0; i <= D; ++i) if(i != low)
                    {
                        vertexSum += simplex[i].first = (simplex[i].first +
                            simplex[low].first) * 0.5;
                        simplex[i].second = f(simplex[i].first);
                    }
                }
            }
        }
    }

    P restartedMinimize(POINT const& initialGuess, int maxIterations = 10000,
        double yPrecision = numeric_limits<double>::epsilon(),
        int maxRepeats = 10, double step = 1)
    {
        P result(initialGuess, numeric_limits<double>::infinity());
        while(maxRepeats--)
        {
            double yOld = result.second;
            result = minimize(result.first, maxIterations, yPrecision, step);
            if(!isELess(result.second, yOld, yPrecision)) break;
        }
        return result;
    }
};

template<typename POINT, typename FUNCTION> struct LineFunction
{
    FUNCTION f;
    POINT x, direction;
    LineFunction(POINT const& theX, POINT const& theDirection,
        FUNCTION const& theF = FUNCTION()): x(theX),
        direction(theDirection), f(theF) {}
    double operator()(double step)const
        {return f(x + direction * step);}
};
template<typename POINT, typename FUNCTION> pair<POINT, double> lineMinimize(
    POINT const& x, POINT const& direction, FUNCTION const& f = FUNCTION(),
    double d = 1)
{
    LineFunction<POINT, FUNCTION> lf(x, direction, f);
    pair<double, double> result =
        minimizeGS(lf, 0, findIntervalBound(lf, 0, d));
    return pair<POINT, double>(x + direction * result.first, result.second);
}

template<typename FUNCTION, typename POINT> pair<POINT, double>
backtrackLineSearch(POINT const& x, POINT const& direction,
    POINT const& gradient, FUNCTION const& f, double y,
    double yPrecision = numeric_limits<double>::epsilon(), double c = 0.0001)
{//dy = step * direction * gradient + O(step * step)`
    double dd = -(direction * gradient), minDecrease = dd * c, yNew;
    if(minDecrease > 0)//ensure decrease
        for(double step = 1; step * (dd + step) >= yPrecision; step /= 2)
        {
            POINT xNew = x + direction * step;
            if(y - (yNew = f(xNew)) >= minDecrease * step)
                return pair<POINT, double>(xNew, yNew);
        }
    return pair<POINT, double>(x, y);
}

template<typename POINT, typename FUNCTION, typename GRADIENT>
    pair<POINT, double> LBFGSMinimize(POINT const& initialGuess,
    FUNCTION const& f, GRADIENT const& g, int maxIterations = 10000,
    double yPrecision = numeric_limits<double>::epsilon(),
    bool useExactSearch = true, int historySize = 8)
{
    Queue<pair<POINT, POINT> > history;
    pair<POINT, double> xy(initialGuess, f(initialGuess));
    POINT grad = g(xy.first), d = -grad;
    while(maxIterations-- > 0)
    {
        double yLast = xy.second;
        pair<POINT, double> xyNew = useExactSearch ?
            lineMinimize(xy.first, d, f) : backtrackLineSearch(xy.first, d,
            grad, f, xy.second, yPrecision);
        if(!isfinite(xyNew.second) ||
            !isELess(xyNew.second, xy.second, yPrecision)) break;
        POINT newGrad = g(xyNew.first);
        if(history.getSize() >= historySize) history.pop();
        history.push(make_pair(xyNew.first - xy.first, newGrad - grad));
        xy = xyNew;
        d = grad = newGrad;//"double recursion" algorithm to update d
        Vector<double> a, p;
        int last = history.getSize() - 1;
        for(int i = last; i >= 0; --i)
        {
            double pi = 1/(history[i].first * history[i].second),
                ai = history[i].first * d * pi;
            d -= history[i].second * ai;
            a.append(ai);
            p.append(pi);
        }//initial Hessian is scaled diagonal
        d *= 1/(history[last].second * history[last].second * p[last]);
        for(int i = 0; i < history.getSize(); ++i)
        {
            double bi = history[i].second * d * p[last - i];
            d += history[i].first * (a[last - i] - bi);
        }
        d *= -1;
    }
    return xy;
}

template<typename FUNCTION> Vector<double> gridMinimize(
    Vector<Vector<double> > const& sets, FUNCTION const& f = FUNCTION())
{
    assert(sets.getSize() > 0);
    long long total = 1;
    for(int i = 0; i < sets.getSize(); ++i)
    {
        assert(sets[i].getSize() > 0);
        total *= sets[i].getSize();
    }
    Vector<double> best;
    double bestScore;
    for(long long i = 0; i < total; ++i)
    {//unrank and eval
        long long rank = i;
        Vector<double> next;
        for(int j = 0; j < sets.getSize(); ++j)
        {
            next.append(sets[j][rank % sets[j].getSize()]);
            rank /= sets[j].getSize();
        }
        double score = f(next);
        if(best.getSize() == 0 || score < bestScore)
        {
            bestScore = score;
            best = next;
        }
    }
    return best;
}

template<typename FUNCTION> Vector<double> randomDiscreteMinimize(
    Vector<Vector<double> > const& sets, FUNCTION const& f = FUNCTION(),
    int evals = 20)
{
    assert(sets.getSize() > 0);
    Vector<double> best;
    double bestScore;
    for(long long i = 0; i < evals; ++i)
    {
        Vector<double> next;
        for(int j = 0; j < sets.getSize(); ++j)
            next.append(sets[j][GlobalRNG.mod(sets[j].getSize())]);
        double score = f(next);
        if(best.getSize() == 0 || score < bestScore)
        {
            bestScore = score;
            best = next;
        }
    }
    return best;
}

template<typename FUNCTION> Vector<double> ILSCompassDiscreteMinimize(
    Vector<Vector<double> > const& sets, FUNCTION const& f = FUNCTION(),
    int remainingEvals = 100)
{
    assert(remainingEvals > 0);
    double bestScore;
    Vector<double> best;
    while(remainingEvals > 0)
    {
        Vector<int> current;
        for(int i = 0; i < sets.getSize(); ++i)
        {
            assert(sets[i].getSize() > 0);
            current.append(sets[i][GlobalRNG.mod(sets[i].getSize())]);
        }
        pair<Vector<double>, pair<double, int> > state =
            compassDiscreteMinimizeHelper(sets, current, f, remainingEvals);
        remainingEvals = state.second.second;
        if(state.second.first < bestScore)
        {
            best = state.first;
            bestScore = state.second.first;
        }
    }
    return best;
}

template<typename FUNCTION> Vector<double> compassDiscreteMinimize(
    Vector<Vector<double> > const& sets, FUNCTION const& f = FUNCTION(),
    int remainingEvals = 100)
{//use median in each set as initial solution
    Vector<int> current;
    for(int i = 0; i < sets.getSize(); ++i)
    {
        assert(sets[i].getSize() > 0);
        current.append(sets[i].getSize()/2);
    }
    return compassDiscreteMinimizeHelper(sets, current, f,
        remainingEvals).first;
}
//assumes set values are in sorted (or reverse sorted) order!
template<typename FUNCTION> pair<Vector<double>, pair<double, int> >
    compassDiscreteMinimizeHelper(Vector<Vector<double> > const& sets,
    Vector<int> current, FUNCTION const& f = FUNCTION(),
    int remainingEvals = 100)
{//start with medians
    Vector<double> best;
    for(int i = 0; i < sets.getSize(); ++i)
    {
        assert(0 <= current[i] && current[i] < sets[i].getSize());
        best.append(sets[i][current[i]]);
    }
    double bestScore = f(best);
    Vector<int> preferredSign(sets.getSize(), 1);
    for(bool isOpt = false, changedSign = false; !isOpt;)
    {
        isOpt = true;
        for(int i = 0; i < sets.getSize(); ++i)
        {
            int next = current[i] + preferredSign[i];
            if(0 <= next && next < sets[i].getSize())
            {
                if(remainingEvals-- < 1)
                    return make_pair(best, make_pair(bestScore, 0));
                best[i] = sets[i][next];
                double score = f(best);
                if(score < bestScore)
                {
                    current[i] = next;
                    bestScore = score;
                    isOpt = false;
                }
                else
                {
                    best[i] = sets[i][current[i]];
                    preferredSign[i] *= -1;
                }
            }
            else preferredSign[i] *= -1;
        }
        if(isOpt){if(!changedSign) isOpt = false;}
        else changedSign = false;
    }
    return make_pair(best, make_pair(bestScore, remainingEvals));
}

double RMRate(int i){return 1/pow(i + 1, 0.501);}

template<typename POINT, typename FUNCTION> POINT SPSA(POINT x,
    FUNCTION const& f, int maxEvals = 10000, double initialStep = 1)
{
    POINT direction = x;
    for(int i = 0, D = x.getSize(); i < maxEvals/2; ++i)
    {
        for(int j = 0; j < D; ++j) direction[j] =
            GlobalRNG.next() % 2 ? 1 : -1;
        double step = initialStep/pow(i + 1, 0.101), temp = RMRate(i) *
            (f(x + direction * step) - f(x - direction * step))/2;
        if(!isfinite(temp)) break;
        for(int j = 0; j < D; ++j) x[j] -= temp/direction[j];
    }
    return x;
}
template<typename POINT, typename FUNCTION> pair<POINT, double> metaSPSA(
    POINT x, FUNCTION const& f, int spsaEvals = 100000, int estimateEvals =
    100, double step = pow(2, 10), double minStep = pow(2, -20))
{
    pair<POINT, double> xy(x, numeric_limits<double>::infinity());
    for(; step > minStep; step /= 2)
    {
        if(isfinite(xy.second)) x = SPSA(xy.first, f, spsaEvals, step);
        double sum = 0;
        for(int i = 0; i < estimateEvals; ++i) sum += f(x);
        if(sum/estimateEvals < xy.second)
        {
            xy.first = x;
            xy.second = sum/estimateEvals;
        }
    }
    return xy;
}

template<typename POINT, typename FUNCTION> pair<POINT, double>
    DSRS(POINT const& initialGuess, FUNCTION const& f,
    double step = 1, double factor = 0.8, int maxFEvals = 10000000,
    double yPrecision = numeric_limits<double>::epsilon())
{
    pair<POINT, double> xy(initialGuess, f(initialGuess));
    for(double dd = 0; --maxFEvals && step * (dd + step) > yPrecision;)
    {
        POINT direction = initialGuess;//ensure non-zero direction
        for(int j = 0; j < direction.getSize(); ++j) direction[j] =
            GlobalRNG.uniform01() * (GlobalRNG.next() % 2 ? 1 : -1);
        direction *= 1/sqrt(direction * direction);
        double yNew = f(xy.first + direction * step);
        if(isELess(yNew, xy.second, yPrecision))
        {
            dd = (xy.second - yNew)/step;
            xy.first += direction * step;
            xy.second = yNew;
            step *= 2;
        }
        else step *= factor;
    }
    return xy;
}

template<typename POINT, typename FUNCTION> struct ILSNelderMead
{
    typedef NelderMead<POINT, FUNCTION> NM;
    typedef typename NM::P P;
    struct Move
    {
        NM &nm;
        P current, best;
        int maxIterations;
        double precision;
        void localSearchBest()
            {current = nm.minimize(current.first, maxIterations, precision);}
        void bigMove()
        {
            for(int i = 0; i < current.first.getSize(); ++i)
                current.first[i] = GlobalRNG.cauchy(0, 1);
        }
        void updateBest()
            {if(best.second > current.second) best = current;}
    };
    NM nelderMead;
    ILSNelderMead(FUNCTION const& theFunction = FUNCTION()):
        nelderMead(theFunction) {}
    P minimize(POINT const& initialGuess, int maxJumps =
        1000, int maxIterations = 1000, double precision = 0.001)
    {
        P initial(initialGuess, numeric_limits<double>::max());
        Move move = {nelderMead, initial, initial, maxIterations, precision};
        iteratedLocalSearch(move, maxJumps);
        return move.best;
    }
};

template<typename POINT> double boxVolume(pair<POINT, POINT>const& box)
{
    double result = 1;
    for(int i = 0; i < box.first.getSize(); ++i)
        result *= box.second[i] - box.first[i];
    return result;
}
template<typename POINT, typename TEST, typename FUNCTION>
double MonteCarloIntegrate(pair<POINT, POINT> const& box, int n,
    TEST const& isInside = TEST(), FUNCTION const& f = FUNCTION())
{
    double sum = 0;
    for(int i = 0; i < n; ++i)
    {
        POINT point;
        for(int j = 0; j < point.getSize(); ++j)
            point[j] = GlobalRNG.uniform(box.first[j], box.second[j]);
        if(isInside(point)) sum += f(point);
    }
    return boxVolume(box) * sum/n;
}

template<typename FUNCTION> class Trapezoid
{
    double sum;
    int nIntervals;
public:
    Trapezoid(): nIntervals(1) {}
    double addLevel(FUNCTION const& f, double xLeft, double xRight)
    {
        double dX = (xRight - xLeft)/nIntervals;
        if(nIntervals == 1) sum = (f(xLeft) + f(xRight))/2;
        else
        {//nIntervals/2 is the number of added points spaced 2*dX apart
            double x = xLeft + dX;
            for(int i = 0; i < nIntervals/2; ++i, x += 2 * dX) sum += f(x);
        }
        nIntervals *= 2;
        return dX * sum;
    }
    static double integrate(FUNCTION const& f, double xLeft, double xRight,
        double maxXError = 0, int maxIterations = 20, int minIterations = 5)
    {
        Trapezoid t;
        double result = t.addLevel(f, xLeft, xRight), oldResult = 0;
        while(--maxIterations > 0 && (--minIterations > 0 ||
            abs(result - oldResult) >= maxXError))
        {
            oldResult = result;
            result = t.addLevel(f, xLeft, xRight);
        }
        return result;
    }
};

class DynamicLinearInterpolation
{
    Treap<double, double> values;
public:
    double findMin()
    {
        assert(!values.isEmpty());
        return values.findMin()->key;
    }
    double findMax()
    {
        assert(!values.isEmpty());
        return values.findMax()->key;
    }
    double evaluate(double x)
    {
        assert(x >= findMin() && x <= findMax());
        double *y = values.find(x);
        if(y) return *y;
        Treap<double, double>::NodeType* left = values.predecessor(x),
            *right = values.successor(x);
        return left->value + (right->value - left->value) *
            (x - left->key)/(right->key - left->key);
    }
    void remove(double x){values.remove(x);}
    bool contains(double x){return values.find(x);}
    void insert(double x, double y){values.insert(x, y);}
};

template<typename ITEM> class Poly
{
    Vector<ITEM> a;
    void trim()
        {while(getSize() > 1 && a.lastItem() == ITEM()) a.removeLast();}
public:
    Poly(Vector<ITEM> const& coefs = Vector<ITEM>(1, ITEM())):
        a(coefs){trim();}
    int getSize()const{return a.getSize();}//degree is size + 1
    void increaseDegree(ITEM const& coef)
    {
        a.append(coef);
        trim();
    }
    ITEM& operator[](int i)
    {
        assert(i >= 0 && i < getSize());
        return a[i];
    }
    ITEM const& operator[](int i)const
    {
        assert(i >= 0 && i < getSize());
        return a[i];
    }
    ITEM operator()(ITEM const& x)const
    {
        ITEM result = ITEM();
        for(int i = getSize() - 1; i >= 0; --i) result = result * x + a[i];
        return result;
    }
    Poly differentiate()const
    {
        Poly result = *this;
        for(int i = 1; i < result.getSize(); ++i)
            result.a[i - 1] = i * result.a[i];
        if(result.getSize() > 0) result.a.removeLast();
        return result;
    }
    Poly integrate()const
    {
        Poly result = *this;
        result.a.append(ITEM());
        for(int i = 0; i < result.getSize() - 1; ++i)
            result.a[i + 1] = result.a[i]/(i + 1);
        return result;
    }
    Poly& operator+=(Poly const& rhs)
    {
        for(int i = 0; i < min(getSize(), rhs.getSize()); ++i) a[i] += rhs[i];
        for(int i = getSize(); i < rhs.getSize(); ++i) a.append(rhs[i]);
        trim();
        return *this;
    }
    Poly operator+(Poly const& rhs)
    {
        Poly result = *this;
        return result += rhs;
    }
    Poly& operator>>=(int k)
    {
        if(k > 0)
        {
            for(int i = k; i < getSize(); ++i) a[i - k] = a[i];
            while(k-- && getSize()) a.removeLast();
            if(getSize() == 0) a.append(ITEM());
        }
        return *this;
    }
    Poly operator>>(int k)
    {
        Poly result = *this;
        return result >>= k;
    }
    Poly& operator<<=(int k)
    {
        if(k > 0)
        {
            for(int i = 0; i < k; ++i) a.append(ITEM());
            for(int i = 0; i < getSize() - k; ++i)
                a[getSize() - 1 - i] = a[getSize() - k - 1 - i];
            for(int i = 0; i < k; ++i) a[i] = ITEM();
        }
        return *this;
    }
    Poly operator<<(int k)
    {
        Poly result = *this;
        return result <<= k;
    }
    Poly& operator*=(ITEM const& rhs)
    {
        a *= rhs;
        return *this;
    }
    Poly operator*(ITEM const& rhs)
    {
        Poly result = *this;
        return result *= rhs;
    }
    Poly& operator-=(Poly const& rhs)
    {
        *this *= -1;
        *this += rhs;
        return *this *= -1;
    }
    Poly operator-(Poly const& rhs)
    {
        Poly result = *this;
        return result -= rhs;
    }
    Poly& operator*=(Poly const& rhs)
    {
        Poly temp = *this <<= rhs.getSize();
        *this *= rhs.a.lastItem();
        for(int i = rhs.getSize() - 2; i >= 0; --i)
        {
            temp >>= 1;
            *this += temp * rhs[i];
        }
        return *this;
    }
    Poly operator*(Poly const& rhs)
    {
        Poly result = *this;
        return result *= rhs;
    }
    void debug()
    {
        DEBUG("POLY BEGIN");
        for(int i = 0; i < getSize(); ++i) DEBUG(a[i]);
        DEBUG("POLY END");
    }
};


/*
use basic dct to find coeffs
then have integration formulas and Clenshaw to evaluate!
assume {-1, 1} and let user do inf transforms?
have getpoly method -> may overflow
error too small if at next ceb point prev poly does well or diff between prev and pre-prev too small?
or coefs -> 0 ? better! on 2^n points only!
*/

struct ChebFunction
{
    Vector<double> ci;
    ChebFunction(Vector<double> const& theCi): ci(theCi)
        {assert(ci.getSize() > 0);}
    template<typename FUNCTION> ChebFunction(FUNCTION const& f, int n)
    {//bar formula need storing both x and y but C formula doesn't
        assert(n > 0);
        Vector<double> fx;
        for(int i = 0; i < n; ++i)
        {
            double x = cos(Random<>::PI() * (i + 0.5)/n);
            fx.append(f(x));
            //DEBUG(x);
            //DEBUG(fx[i]);
        }
        for(int j = 0; j < n; ++j)
        {
            double c = 0;
            for(int i = 0; i < n; ++i)
            {
                double x = j * Random<>::PI() * (i + 0.5)/n;
                c += fx[i] * cos(x);
            }
            ci.append(c * (j > 0 ? 2 : 1)/n);
            //DEBUG(ci[j]);
        }
    }
    double operator()(double x)
    {
        assert(x >= -1 && x <= 1);
        double d = 0, dd = 0, y2 = 2 * x;
        for(int i = ci.getSize() - 1; i > 0; --i)
        {
            double sv = d;
            d = y2 * d - dd + ci[i];
            dd = sv;
        }
        return x * d - dd + ci[0];
    }
    ChebFunction integral(double FM1 = 0)
    {
        Vector<double> result;
        result.append(FM1);
        for(int i = 1; i - 1 < ci.getSize(); ++i)
            result.append(((i - 1 > 0 ? 1 : 2) * ci[i - 1] -
                (i + 1 > ci.getSize() - 1 ? 0 : ci[i + 1]))/2/i);
        ChebFunction cf(result);
        cf.ci[0] -= cf(-1);
        return cf;
    }
    ChebFunction derivative()
    {
        int n = ci.getSize() - 1;
        Vector<double> result(n, 0);
        for(int i = n; i > 0; --i) result[i - 1] =
            (i + 1 > n - 1 ? 0 : result[i + 1]) + 2 * i * ci[i];
        result[0] /= 2;
        return ChebFunction(result);
    }
    double integrateM11()
    {
        double result = 0;
        for(int i = 0; i < ci.getSize(); i += 2)
            result += 2 * ci[i]/(1 - i * i);
        return result;
    }
    Poly<double> poly()
    {
        Poly<double> result;
        Vector<Poly<double> > polys;
        polys.append(Poly<double>(Vector<double>(1, 1)));
        Vector<double> temp(2, 0);
        temp[1] = 1;
        polys.append(Poly<double>(temp));
        for(int i = 0; i < ci.getSize(); ++i)
        {
            if(i > 1) polys.append(
                ((polys[i - 1] * 2) << 1) - polys[i - 2]);
            result += polys[i] * ci[i];
        }
        return result;
    }
};

template<typename FUNCTION> class ScaledFunctionM11
{
    FUNCTION f;
    double a, b;
public:
    ScaledFunctionM11(double theA, double theB, FUNCTION const& theF =
        FUNCTION()): f(theF), a(theA), b(theB) {assert(a < b);}
    double operator()(double u)const{return f(((b - a) * u + a + b)/2);}
};

template<typename FUNCTION> class ScaledFunctionAB
{
    FUNCTION f;
    double a, b;
public:
    ScaledFunctionAB(double theA, double theB, FUNCTION const& theF =
        FUNCTION()):f(theF), a(theA), b(theB) {assert(a < b);}
    double operator()(double x)const{return (2 * x + a + b)/(b - a);}
};

template<typename FUNCTION> double integrate(FUNCTION const& f, double a,
    double b, double maxXError = numeric_limits<double>::epsilon(),
    int maxEvals = 1025, int minIterations = 32)
{
    int n = minIterations;
    ChebFunction che(ScaledFunctionM11<FUNCTION>(a, b, f), n);
    double result = che.integrateM11();
    while(n < maxEvals)
    {
        ChebFunction che2(ScaledFunctionM11<FUNCTION>(a, b, f), n *= 2);
        double result2 = che2.integrateM11();
        if(isEEqual(result, result2, maxXError)) return result2;
        result = result2;
        che = che2;
    }
    return result;
}

typedef complex<double> Complex;
Vector<Complex> fastFourierTransform(Vector<Complex> const& x)
{//FFT only makes sense for powers of 2
    int n = x.getSize(), b = lgFloor(n);
    assert(isPowerOfTwo(n));
    Vector<Complex> result(n);
    for(int i = 0; i < n; ++i) result[reverseBits(i, b)] = x[i];
    for(int s = 1; s <= b; ++s)
    {
        int m = twoPower(s);
        Complex wm = exp((2 * Random<>::PI()/m) * Complex(0, 1));
        for(int k = 0; k < n; k += m)
        {
            Complex w(1, 0);
            for(int j = 0; k < m/2; ++j)
            {
                Complex t = w * result[k + j + m/2];
                Complex u = result[k + j];
                result[k + j] = u + t;
                result[k + j + m/2] = u - t;
                w *= wm;
            }
        }
    }
    return result;
}

}
#endif
