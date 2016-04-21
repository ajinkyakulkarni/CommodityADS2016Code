#include "ChainingHashTable.h"
#include "LinearProbingHashTable.h"
#include "BloomFilter.h"
#include "../RandomNumberGeneration/Statistics.h"
#include <iostream>
#include <cmath>
using namespace igmdk;

struct Fat2
{
	enum{SIZE = 10};
	int array[SIZE];
	Fat2(int last)
	{
		for(int i = 1; i < SIZE; ++i)
		{
			array[i] = i;
		}
		array[0] = last;
	}
	bool operator==(Fat2 const& rhs)const
	{
		for(int i = 0; i < SIZE; ++i)
		{
			if(array[i] != rhs.array[i]) return false;
		}
		return true;
	}
	bool operator<(Fat2 const& rhs)const
	{
		for(int i = 0; i < SIZE; ++i)
		{
			if(array[i] < rhs.array[i]) return true;
			if(array[i] > rhs.array[i]) return false;
		}
		return false;
	}
	int getSize()const{return SIZE;}
	int const operator[](int i)const{return array[i];}
};

struct Action
{
    void operator()(int const& key, int const& item)
    {
        cout << key << endl;
    }
};

void timeRT()
{
	ChainingHashTable<Fat2, int, BHash<FNVHash> > t;
	//LinearProbingHashTable<int, int, BHash<FNVHash> >t;
	int N = 150000;
	for(int i = 0; i < N; ++i)
	{
		t.insert(i,i);
	}
	for(int j = 0; j < 100; ++j)
	{
		for(int i = 0; i < N; ++i)
		{
			assert(t.find(i));
			assert(*t.find(i) == i);
			//t.remove(i);
		}
	}

}

void timeRT2()
{
    //PrimeHash h;
    //TableHash h;
    //XorshiftHash64 h;
    //XorshiftHash32 h;
    //MHash<XorshiftHash64> h(3573489593);
    //MHash<XorshiftHash32> h(3573489593);
    //MHash<TableHash> h(3573489593);
    MHash<FNVHash> h(3573489593u);
    //BHash h(6);
    //BHash2<XorshiftHash64> h(31);
    //BHash2<XorshiftHash32> h(31);
    //BHash2<TableHash> h(31);
    //BHash2<PrimeHash> h(31);
    //FairHash h;
    //EHash32<MHash<PrimeHash> > h(3573489593);
    //BHash2<PrimeHash> h;
    //EHash32<BHash<FairHash> > h;
    unsigned int sum = 0;
    for(int i = 0; i < 1500000000; ++i)
	{
		sum += h.hash(i);
	}
	DEBUG(sum);
}

struct FunctionTester
{
    void operator()()const
    {
        timeRT();
    }
};

void DDDChaining()
{
    ChainingHashTable<int, int> chainingH0to9;
    for(int i = 0; i < 10; ++i)
	{
		chainingH0to9.insert(i, i);
	}
    cout << "breakpoint" << endl;
}

void DDDLinearProbing()
{
    LinearProbingHashTable<int, int> linearProbingH0to9;
    for(int i = 0; i < 10; ++i)
	{
		linearProbingH0to9.insert(i, i);
	}
    cout << "breakpoint" << endl;
}
/*
void DDDBloomFilter()
{
    BloomFilter<int, int> linearProbingH0to9;
    for(int i = 0; i < 10; ++i)
	{
		linearProbingH0to9.insert(i, i);
	}
}*/

int main()
{
    DDDChaining();
    DDDLinearProbing();
    timeRT();
	/*NormalSummary result = MonteCarlo::simulate(SpeedTester<FunctionTester>(), 1);
    DEBUG(result.minimum);
    DEBUG(result.maximum);
    DEBUG(result.mean);
    DEBUG(result.variance);
    DEBUG(result.confidence9973());*/
	return 0;
}
