#include "UnionFind.h"
#include <iostream>
using namespace igmdk;

int main()
{
	UnionFind uf(10);
	cout << uf.areEquivalent(4, 8) << endl;
	cout << uf.areEquivalent(4, 9) << endl;
	uf.join(4,8);
	cout << uf.areEquivalent(4, 8) << endl;
	cout << uf.areEquivalent(4, 9) << endl;
	uf.increaseSize(20);
	cout << uf.areEquivalent(4, 8) << endl;
	cout << uf.areEquivalent(4, 9) << endl;

	IntervalSetUnion iu;

	iu.split(5);
	DEBUG(iu.find(2));
	iu.split(3);
	DEBUG(iu.find(2));
	iu.merge(3);
	DEBUG(iu.find(2));
	return 0;
}
