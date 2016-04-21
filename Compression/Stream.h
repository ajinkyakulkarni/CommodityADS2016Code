#ifndef STREAM_H
#define STREAM_H
#include "../Utils/Bitset.h"
#include "../Utils/Vector.h"
namespace igmdk{

Vector<unsigned char> ReinterpretEncode(unsigned long long n, int size)
{
    Vector<unsigned char> result;
    while(size-- > 0)
    {
        result.append(n % 256);
        n /= 256;
    }
    return result;
}
unsigned long long ReinterpretDecode(Vector<unsigned char> const& code)
{
    unsigned long long n = 0, base = 1;
    for(int i = 0; i < code.getSize(); ++i)
    {
        n += base * (code[i] % 256);
        base *= 256;
    }
    return n;
}

struct Stream
{
    unsigned long long position;
    Stream(): position(0) {}
};
struct BitStream : public Stream
{
    Bitset<unsigned char> bitset;//unsigned char for portability
    enum{B = 8};
    BitStream() {}
    BitStream(Bitset<unsigned char> const& aBitset): bitset(aBitset) {}
    BitStream(Vector<unsigned char> const& vector): bitset(vector) {}
    void writeBit(bool value){bitset.append(value);}
    bool readBit()
    {
        assert(bitsLeft());
        return bitset[position++];
    }
    void writeByte(unsigned char byte){writeValue(byte, B);}
    unsigned char readByte(){return readValue(B);}
    void output()const{bitset.output();}
    void writeValue(unsigned long long value, int bits)
        {bitset.appendValue(value, bits, true);}
    unsigned long long readValue(int bits)
    {
        assert(bits <= bitsLeft());
        position += bits;
        return bitset.getBitReversedValue(position - bits, bits);
    }
    unsigned long long bitsLeft()const{return bitset.getSize() - position;}
    unsigned long long bytesLeft()const{return bitsLeft()/B;}
};
}
#endif
