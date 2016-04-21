#include "Bitset.h"
using namespace igmdk;

void DDDBitset()
{
    Bitset<unsigned char> BitsetChar19Every4(19);
	for(int i = 0; i < 19; i += 4)
	{
		BitsetChar19Every4.set(i, true);
	}
	cout << "breakpoint" << endl;
}

void DDDNBitVector()
{
    KBitVector<4, unsigned char> Vector4BitChar8to12;
	for(int i = 0; i < 5; i += 1)
	{
		Vector4BitChar8to12.append(i + 8);
	}
	cout << "breakpoint" << endl;
}

int main()
{
    DDDBitset();
    DDDNBitVector();
    for(int i = 0; i < 5; ++i) DEBUG(nextPowerOfTwo(i));

	Bitset<unsigned char> bs(21);
	bs.output();
	for(int i = 0; i < 21; i+=3)
	{
		bs.set(i, true);
		bs.output();
	}
	bs.set(3, false);
	bs.output();
	bs <<= 9;
	bs.output();
	bs >>= 9;
	bs.output();
	bs.flip();
	bs.output();

    KBitVector<5> x;
    int N = 100000;
    for(int i = 0; i < N; ++i) x.append(i);
    for(int i = 0; i < N; ++i) x.set(i, i);
    for(int i = 0; i < N; ++i) assert(x[i] == i % 32);

    DEBUG(reverseBits<unsigned int>(7, 3));
	return 0;
}
