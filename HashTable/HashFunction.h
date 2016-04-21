#ifndef HASH_FUNCTION_H
#define HASH_FUNCTION_H
#include <string>
#include "../RandomNumberGeneration/Random.h"
#include "../Utils/Bits.h"
namespace igmdk{

struct FairHash
{
    unsigned int hash(unsigned int x){return x;}const
    unsigned int hash(unsigned char* array, int size)const
    {
        unsigned int result = 0;
        for(int i = 0; i < min(size, 4); ++i)
            result = (result << 8) | array[i];
        return result;
    }
    template<typename NUMBER> unsigned int hash(NUMBER* array, int size)const
        {return hash((unsigned char*)array, size * sizeof(NUMBER));}
};

struct LCGHash
{
    unsigned int hash(unsigned int x)const{return 1099087573u * x;}
    template<typename NUMBER> unsigned int hash(NUMBER* array, int size)const
    {
        unsigned int sum = 0;
        for(int i = 0; i < size; ++i) sum = 1099087573u * (sum + array[i]);
        return sum;
    }
};

struct FNVHash
{
    template<typename POD> unsigned int hash(POD const& x)const
        {return hash((unsigned char*)&x, sizeof(x));}
    unsigned int hash(unsigned char* array, int size)const
    {
        unsigned int sum = 2166136261u;
        for(int i = 0; i < size; ++i) sum = (sum * 16777619) ^ array[i];
        return sum;
    }
    template<typename NUMBER> unsigned int hash(NUMBER* array, int size)
        {return hash((unsigned char*)array, size * sizeof(NUMBER));}
};

struct FNVHash64
{
    template<typename POD> unsigned long long hash(POD const& x)const
        {return hash((unsigned char*)&x, sizeof(x));}
    unsigned long long hash(unsigned char* array, int size)const
    {
        unsigned long long sum = 14695981039346656037ull;
        for(int i = 0; i < size; ++i)
            sum = (sum * 1099511628211ull) ^ array[i];
        return sum;
    }
    template<typename NUMBER>unsigned long long hash(NUMBER* array, int size)
        const{return hash((unsigned char*)array, size * sizeof(NUMBER));}
};

class Xorshift32Hash
{
    unsigned int seed;
public:
    Xorshift32Hash(): seed(GlobalRNG.next()) {}
    unsigned int hash(unsigned int x)const
        {return Xorshift::transform(seed + x);}
    template<typename NUMBER> unsigned int hash(NUMBER* array, int size)const
    {
        unsigned int sum = seed;
        for(int i = 0; i < size; ++i)
            sum = Xorshift::transform(sum + array[i]);
        return sum;
    }
};

class Xorshift64Hash
{
    unsigned long long seed;
public:
    Xorshift64Hash(): seed(GlobalRNG.next()) {}
    unsigned int hash(unsigned long long x)const
        {return QualityXorshift64::transform(seed + x);}
    template<typename NUMBER> unsigned int hash(NUMBER* array, int size)const
    {
        unsigned long long sum = seed;
        for(int i = 0; i < size; ++i)
            sum = QualityXorshift64::transform(sum + array[i]);
        return sum;
    }
};

class PrimeHash
{
    static unsigned int const PRIME = (1ull << 32) - 5;
    unsigned long long seed;//could be > PRIME but that's ok
public:
    PrimeHash(): seed(GlobalRNG.next()) {}
    unsigned int hash(unsigned int x)const{return seed * x % PRIME;}
    template<typename NUMBER> unsigned int hash(NUMBER* array, int size)const
    {//numbers in the array must fit into unsigned int
        unsigned long long sum = 0, a = seed;
        for(int i = 0; i < size; ++i)
        {
            sum += a * array[i];
            a = Xorshift::transform(a);
        }//possible overflow but that's ok
        return sum % PRIME;
    }
};

class PrimeHash2
{
    static unsigned int const PRIME = (1ull << 32) - 5;
    unsigned long long seed;//could be > PRIME but that's ok
public:
    PrimeHash2(): seed(GlobalRNG.next() % PRIME) {}
    unsigned int hash(unsigned int x)const{return seed * x % PRIME;}
    template<typename NUMBER> unsigned int hash(NUMBER* array, int size)const
    {//numbers in the array must fit into unsigned int
        unsigned long long sum = 0;
        for(int i = 0; i < size; ++i) sum = (sum + seed * array[i]) % PRIME;
        return sum;
    }
};

template<typename HASHER = PrimeHash> struct EHash32
{
    HASHER h;
    EHash32() {}
    EHash32(unsigned long long data): h(data) {}
    typedef unsigned int WORD;
    WORD hash(WORD x)const{return h.hash(x);}
    WORD hash(int x)const{return h.hash((WORD)x);}
    WORD hash(short x)const{return h.hash((WORD)x);}
    WORD hash(unsigned short x)const{return h.hash((WORD)x);}
    template<typename POD> WORD hash(POD x)const
    {
        int rem = sizeof(x) % sizeof(WORD);
        if(rem == 0) return hash((WORD*)&x, sizeof(x)/sizeof(WORD));
        else if(rem == 2) return hash((unsigned short*)&x,
            sizeof(x)/sizeof(unsigned short));
        else return hash((unsigned char*)&x, sizeof(x));
    }
    WORD hash(float x)const
    {
        union{float a; WORD b;} c;
        c.a = x;
        return h.hash(c.b);
    }
    template<typename NUMBER> unsigned int hash(NUMBER* array, int size)const
        {return h.hash(array, size);}
};

template<typename HASHER = PrimeHash> struct EHash64
{
    HASHER h;
    typedef unsigned long long WORD;
    WORD hash(WORD x)const{return h.hash(x);}
    WORD hash(long long x)const{return h.hash((WORD)x);}
    WORD hash(int unsigned x)const{return h.hash((WORD)x);}
    WORD hash(int x)const{return h.hash((WORD)x);}
    WORD hash(short x)const{return h.hash((WORD)x);}
    WORD hash(unsigned short x)const{return h.hash((WORD)x);}
    template<typename POD> WORD hash(POD x)const
    {
        int rem = sizeof(x) % sizeof(WORD);
        if(rem == 0) return hash((WORD*)&x, sizeof(x)/sizeof(WORD));
        else if(rem == 4)
            return hash((unsigned int*)&x, sizeof(x)/sizeof(unsigned int));
        else if(rem == 2) return hash((unsigned short*)&x,
            sizeof(x)/sizeof(unsigned short));
        else return hash((unsigned char*)&x, sizeof(x));
    }
    unsigned int hash(float x)const
    {
        union{float a; WORD b;}c;
        c.a = x;
        return h.hash(c.b);
    }
    template<typename NUMBER> unsigned int hash(NUMBER* array, int size)const
        {return h.hash(array, size);}
};

class TableHash
{
    unsigned int table[256];
public:
    TableHash(){for(int i = 0; i < 256; ++i) table[i] = GlobalRNG.next();}
    unsigned int hash(unsigned char* array, int size)const
    {
        unsigned int result = 0;
        for(int i = 0; i < size; ++i) result ^= table[array[i]];
        return result;
    }
    template<typename NUMBER>unsigned long long hash(NUMBER* array, int size)
        const{return hash((unsigned char*)array, size * sizeof(NUMBER));}
    template<typename POD> unsigned int hash(POD const& x)
        const{return hash((unsigned char*)&x, sizeof(x));}
    unsigned int update(unsigned int currentHash, unsigned char byte)
        const{return currentHash ^ table[byte];}//for both add and remove
};

template<typename HASHER = PrimeHash> class MHash
{
    unsigned int m;
    HASHER h;
public:
    MHash(unsigned int theM = 1ull << 31): m(theM) {}
    template<typename POD> unsigned int hash(POD const& x)const
        {return h.hash(x) % m;}
    template<typename NUMBER> unsigned int hash(NUMBER* array, int size)const
        {return h.hash(array, size) % m;}
};
template<typename HASHER = PrimeHash> class BHash
{
    unsigned int mask;
    HASHER h;
public:
    BHash(unsigned int b = 31): mask(Bits::lowerMask(b)) {}
    template<typename POD> unsigned int hash(POD const& x)const
        {return h.hash(x) & mask;}
    template<typename POD> unsigned int hash(POD* array, int size)const
        {return h.hash(array, size) & mask;}
};

class BUHash
{
    unsigned int a, wLB;
    BHash<PrimeHash> h;
public:
    BUHash(unsigned int b = 31): a(GlobalRNG.next() | 1),
        wLB(numeric_limits<unsigned int>::digits - b), h(b) {}
    unsigned int hash(unsigned int x)const{return (a * x) >> wLB;}
    template<typename POD> unsigned int hash(POD* array, int size)const
        {return h.hash(array, size);}
};

template<typename HASHER = EHash32<BUHash> > struct DataHash
{
    HASHER h;
    unsigned int hash(string const& item)const
        {return h.hash(item.c_str(), item.size());}
    template<typename VECTOR> unsigned int hash(VECTOR const& item)const
        {return h.hash(item.getArray(), item.getSize());}
};

}
#endif
