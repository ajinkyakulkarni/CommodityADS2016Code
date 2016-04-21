#ifndef BITSET_H
#define BITSET_H
#include "../Utils/Vector.h"
#include "../Utils/Utils.h"
#include "Bits.h"
using namespace std;
namespace igmdk{

template<typename WORD = unsigned long long> class Bitset
{
    enum{B = numeric_limits<WORD>::digits, SHIFT = B - 1};
    unsigned long long bitSize;
    mutable Vector<WORD> storage;
    void zeroOutRemainder()const
        {storage.lastItem() &= Bits::upperMask(B - lastWordBits());}
    bool get(int i)const
    {
        if(!(i >= 0 && i < bitSize)) DEBUG(i);
        assert(i >= 0 && i < bitSize);
        return Bits::get(storage[i/B], SHIFT - i % B);
    }
    unsigned long long wordsNeeded()const{return ceiling(bitSize, B);}
public:
    Bitset(unsigned long long initialSize = 0):
       bitSize(initialSize), storage(wordsNeeded(), 0){}
    Bitset(Vector<WORD> const& vector):
        storage(vector), bitSize(B * vector.getSize()){}
    int lastWordBits()const
    {
        assert(bitSize > 0);
        int result = bitSize % B;
        if(result == 0) result = B;
        return result;
    }
    int garbageBits()const{return B - lastWordBits();}
    Vector<WORD>& getStorage(){return storage;}
    unsigned long long getSize()const{return bitSize;}
    unsigned long long wordSize()const{return storage.getSize();}
    bool operator[](int i)const{return get(i);}
    void append(bool value)
    {
        ++bitSize;
        if(wordSize() < wordsNeeded()) storage.append(0);
        set(bitSize - 1, value);
    }
    void set(int i, bool value = true)
    {
        if(!(i >= 0 && i < bitSize)) DEBUG(i);
        assert(i >= 0 && i < bitSize);
        Bits::set(storage[i/B], SHIFT - i % B, value);
    }
    void removeLast()
    {
        assert(bitSize > 0);
        if(lastWordBits() == 1) storage.removeLast();
        --bitSize;
    }
    void setAll(bool value = true)
    {
        for(int i = 0; i < wordSize(); ++i)
            storage[i] = value ? Bits::FULL : Bits::ZERO;
    }
    bool isZero()
    {
        zeroOutRemainder();
        for(int i = 0; i < wordSize(); ++i) if(storage[i])return false;
        return true;
    }
    unsigned long long getValue(int i, int n)const
    {
        assert(n <= numeric_limits<unsigned long long>::digits && i >= 0 &&
            i + n <= bitSize);
        unsigned long long result = 0;
        for(int j = 0; j < n; ++j) Bits::set(result, n - 1 - j, get(i + j));
        return result;
    }
    unsigned long long getBitReversedValue(int i, int n)const
    {
        assert(n <= numeric_limits<unsigned long long>::digits && i >= 0 &&
            i + n <= bitSize && n > 0);
        i += n - 1;
        int index = i/B, right = SHIFT - i % B;
        unsigned long long result = Bits::getValue(storage[index], right, n);
        for(int m = B - right; m < n; m += B)
            result |= Bits::getValue(storage[--index], 0, n - m) << m;
        return result;
    }
    void setValue(unsigned long long value, int i, int n)
    {
        assert(n <= numeric_limits<unsigned long long>::digits && i >= 0 &&
            i + n <= bitSize && n > 0);
        for(int j = 0; j < n; ++j) set(i + j, Bits::get(value, n - 1 - j));
    }
    void setBitReversedValue(unsigned long long value, int i, int n)
    {
        assert(n <= numeric_limits<unsigned long long>::digits && i >= 0 &&
            i + n <= bitSize && n > 0);
        i += n - 1;
        int index = i/B, right = SHIFT - i % B;
        assert(index >= 0 && index < storage.getSize());
        Bits::setValue(storage[index], value, right, n);
        for(int m = B - right; m < n; m += B)
        {
            --index;
            assert(index >= 0);
            Bits::setValue(storage[index], value >> m, 0, n - m);
        }
    }
    void appendValue(unsigned long long value, int n,
        bool reverseBits = false)
    {
        int start = bitSize;
        bitSize += n;
        int k = wordsNeeded() - wordSize();
        for(int i = 0; i < k; ++i) storage.append(0);
        if(reverseBits) setBitReversedValue(value, start, n);
        else setValue(value, start, n);
    }
    void appendBitset(Bitset const& rhs)
    {
        for(int i = 0; i < rhs.wordSize(); ++i)
            appendValue(rhs.storage[i], B, true);
        bitSize -= B - rhs.lastWordBits();
        if(wordSize() > wordsNeeded()) storage.removeLast();
    }

    bool operator==(Bitset const& rhs)const{return storage == rhs.storage;}
    Bitset& operator&=(Bitset const& rhs)
    {
        for(int i = 0; i < min(wordSize(), rhs.wordSize()); ++i)
            storage[i] &= rhs.storage[i];
        return *this;
    }
    void flip()
        {for(int i = 0; i < wordSize(); ++i) storage[i] = ~storage[i];}

    Bitset& operator|=(Bitset const& rhs)
    {
        for(int i = 0; i < min(wordSize(), rhs.wordSize()); ++i)
            storage[i] |= rhs.storage[i];
        return *this;
    }
    Bitset& operator^=(Bitset const& rhs)
    {
        for(int i = 0; i < min(wordSize(), rhs.wordSize()); ++i)
            storage[i] ^= rhs.storage[i];
        return *this;
    }
    Bitset& operator>>=(int shift)
    {
        if(shift < 0) return operator<<=(-shift);
        int normalShift = shift % bitSize, wordShift = normalShift/B,
            bitShift = normalShift % B;
        if(wordShift > 0)//shift words
            for(int i = 0; i + wordShift < wordSize(); ++i)
            {
                storage[i] = storage[i + wordShift];
                storage[i + wordShift] = 0;
            }
        if(bitShift > 0)//shift bits
        {//little endian shift 00000101|00111000 >>= 4 -> 01010011|10000000
            WORD carry = 0;
            zeroOutRemainder();
            for(int i = wordSize() - 1 - wordShift; i >= 0; --i)
            {
                WORD tempCarry = storage[i] >> (SHIFT - bitShift);
                storage[i] <<= bitShift;
                storage[i] |= carry;
                carry = tempCarry;
            }
        }
        return *this;
    }
    Bitset& operator<<=(int shift)
    {
        if(shift < 0) return operator>>=(-shift);
        int normalShift = shift % bitSize, wordShift = normalShift/B,
            bitShift = normalShift % B;
        if(wordShift > 0)//shift words
            for(int i = wordSize() - 1; i - wordShift >= 0; --i)
            {
                storage[i] = storage[i - wordShift];
                storage[i - wordShift] = 0;
            }
        if(bitShift > 0)//shift bits
        {////little endian shift 01010011|10000000 <<= 4 -> 00000101|00111000
            WORD carry = 0;
            for(int i = wordShift; i < wordSize(); ++i)
            {
                WORD tempCarry = storage[i] << (SHIFT - bitShift);
                storage[i] >>= bitShift;
                storage[i] |= carry;
                carry = tempCarry;
            }
        }
        return *this;
    }
    void output()const
    {
        for(int i = 0; i < bitSize; ++i) cout << get(i);
        cout << endl;
    }
    void reverse()
    {
        for(int i = 0; i < bitSize/2; ++i)
        {
            bool temp = get(i);
            set(i, get(bitSize - 1 - i));
            set(bitSize - 1 - i, temp);
        }
    }
    int popCount()const
    {
        zeroOutRemainder();
        int sum = 0;
        for(int i = 0; i < wordSize(); ++i) sum += popCountWord(storage[i]);
        return sum;
    }
};

template<int N, typename WORD = unsigned long long> class KBitVector
{
    Bitset<WORD> bitset;
public:
    WORD operator[](unsigned long long i){return bitset.getValue(i * N, N);}
    void set(WORD value, unsigned long long i)
        {bitset.setValue(value, i * N, N);}
    void append(WORD value){bitset.appendValue(value, N);}
    unsigned long long getSize(){return bitset.getSize()/N;}
};

}
#endif
