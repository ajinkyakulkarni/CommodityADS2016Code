#include "Misc.h"
#include <cassert>
#include <iostream>
using namespace std;
using namespace igmdk;

void testLRU()
{
    LRUCache<int> c(5);
    for(int i = 0; i < 10; ++i) c.write(i, i);
    for(LRUCache<int>::Iterator i = c.begin(); i != c.end(); ++i)
        DEBUG(i->value);
    for(int i = 5; i < 10; ++i)
    {
        DEBUG(i);
        assert(c.read(i));
        DEBUG(*c.read(i));
        assert(*c.read(i) == i);
    }
}

void DDDLRU()
{
    LRUCache<int> LRU4_0to9(4);
    for(int i = 0; i < 10; ++i) LRU4_0to9.write(i, i);

    cout << "breakpoint" << endl;
}


int main()
{
    DDDLRU();
    testLRU();
    return 0;

    //Combinator perm(2, 4);
    Permutator perm(4);
    for(int j = 0; j < 7; ++j)
    {
        cout << "P" << endl;
        for(int i = 0; i < perm.p.getSize(); ++i) cout << perm.p[i] << " ";
        cout << endl;
        if(perm.next()) break;
    }
    DEBUG("done1");
    perm.advance(0);
    perm.advance(1);
    for(;;)
    {
        cout << "P" << endl;
        for(int i = 0; i < perm.p.getSize(); ++i) cout << perm.p[i] << " ";
        cout << endl;
        if(perm.next()) break;
    }


    DEBUG("done2");

    Combinator comb(4, 6);

    for(int j = 0; j < 7; ++j)
    {
        cout << "C" << endl;
        for(int i = 0; i < comb.c.getSize(); ++i) cout << comb.c[i] << " ";
        cout << endl;
        if(comb.next()) break;
    }
    DEBUG("done3");
    comb.skipAfter(1);
    for(;;)
    {
        cout << "C" << endl;
        for(int i = 0; i < comb.c.getSize(); ++i) cout << comb.c[i] << " ";
        cout << endl;
        if(comb.next()) break;
    }
    DEBUG("done4");
    Partitioner part(4);

    for(;;)
    {
        cout << "Pa" << endl;
        for(int i = 0; i < part.p.getSize(); ++i) cout << part.p[i] << " ";
        cout << endl;
        if(part.next()) break;
    }


	return 0;
}
