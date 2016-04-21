#ifndef POINT_H
#define POINT_H
#include "../Utils/Utils.h"
#include <cmath>
#include <cstdlib>//int version of abs
namespace igmdk{

template<typename KEY, int D = 2> class Point
{
    KEY x[D];
public:
    static int const d = D;
    KEY& operator[](int i){assert(i >= 0 && i < D); return x[i];}
    KEY const& operator[](int i)const{assert(i >= 0 && i < D); return x[i];}
    int getSize()const{return D;}
    Point(){for(int i = 0; i < D; ++i) x[i] = 0;}
    Point(KEY const& x0, KEY const& x1)
    {
        assert(D > 1);
        x[0] = x0;
        x[1] = x1;
    }
    bool operator==(Point const& rhs)const
    {
        for(int i = 0; i < D; ++i) if(x[i] != rhs.x[i]) return false;
        return true;
    }
    Point& operator+=(Point const& rhs)
    {
        for(int i = 0; i < D; ++i) x[i] += rhs.x[i];
        return *this;
    }
    friend Point operator+(Point const& lhs, Point const& rhs)
    {
        Point result = lhs;
        return result += rhs;
    }
    Point& operator*=(double scalar)
    {
        for(int i = 0; i < D; ++i) x[i] *= scalar;
        return *this;
    }
    friend Point operator*(Point const& point, double scalar)
    {
        Point result = point;
        return result *= scalar;
    }
    //new operators
    Point& operator-=(Point const& rhs){return *this += rhs * -1;}
    friend Point operator-(Point const& lhs, Point const& rhs)
    {
        Point result = lhs;
        return result -= rhs;
    }
    Point operator-(){return *this * -1;}
    double operator*(Point const& rhs)const
    {
        double dotProduct = 0;
        for(int i = 0; i < D; ++i) dotProduct += x[i] * rhs[i];
        return dotProduct;
    }
};
typedef Point<double> Point2;

template<typename VECTOR> struct EuclideanDistance
{
    static double iDistanceIncremental(VECTOR const& lhs, VECTOR const& rhs,
        int i)
    {
        double x = lhs[i] - rhs[i];
        return x * x;
    }
    static double distanceIncremental(VECTOR const& lhs, VECTOR const& rhs,
        double bound = numeric_limits<double>::max())
    {
        assert(lhs.getSize() == rhs.getSize());
        double sum = 0;
        for(int i = 0; i < lhs.getSize() && sum < bound; ++i)
            sum += iDistanceIncremental(lhs, rhs, i);
        return sum;
    }
    struct Distance
    {
        double operator()(VECTOR const& lhs, VECTOR const& rhs)const
            {return sqrt(distanceIncremental(lhs, rhs));}
    };
    struct DistanceIncremental
    {
        double operator()(VECTOR const& lhs, VECTOR const& rhs)const
            {return distanceIncremental(lhs, rhs);}
        double operator()(VECTOR const& lhs, VECTOR const& rhs, int i)const
            {return iDistanceIncremental(lhs, rhs, i);}
        double operator()(double bound, VECTOR const& lhs, VECTOR const& rhs)
            const{return distanceIncremental(lhs, rhs, bound);}
    };
};

template<typename VECTOR> struct MaxDistance
{
    static double iDistanceIncremental(
        VECTOR const& lhs, VECTOR const& rhs, int i)
        {return abs(lhs[i] - rhs[i]);}
    static double distanceIncremental(VECTOR const& lhs,
        VECTOR const& rhs, double bound = numeric_limits<double>::max())
    {
        assert(lhs.getSize() == rhs.getSize());
        double sum = 0;
        for(int i = 0; i < lhs.getSize() && sum < bound; ++i)
            sum = max(sum, iDistanceIncremental(lhs, rhs, i));
        return sum;
    }
    struct Distance
    {
        double operator()(VECTOR const& lhs, VECTOR const& rhs)const
            {return distanceIncremental(lhs, rhs);}
    };
    struct DistanceIncremental
    {
        double operator()(VECTOR const& lhs, VECTOR const& rhs)const
            {return distanceIncremental(lhs, rhs);}
        double operator()(VECTOR const& lhs, VECTOR const& rhs, int i)const
            {return iDistanceIncremental(lhs, rhs, i);}
        double operator()(double bound, VECTOR const& lhs, VECTOR const& rhs)
            const{return distanceIncremental(lhs, rhs, bound);}
    };
};

}//end namespace
#endif
