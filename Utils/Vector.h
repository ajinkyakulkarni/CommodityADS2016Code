#ifndef VECTOR_H
#define VECTOR_H
#include "Utils.h"
namespace igmdk{

template<typename ITEM> class Vector
{
    enum{MIN_CAPACITY = 8};
    int capacity, size;
    ITEM* items;
    void resize()
    {
        ITEM* oldItems = items;
        capacity = max(2 * size, int(MIN_CAPACITY));
        items = rawMemory<ITEM>(capacity);
        for(int i = 0; i < size; ++i) new(&items[i])ITEM(oldItems[i]);
        rawDestruct(oldItems, size);
    }
public:
    ITEM* getArray(){return items;}
    ITEM* const getArray()const{return items;}
    int getSize()const{return size;}
    ITEM& operator[](int i)
    {
        assert(i >= 0 && i < size);
        return items[i];
    }
    ITEM const& operator[](int i)const
    {
        assert(i >= 0 && i < size);
        return items[i];
    }

    Vector(ITEM* const array = 0, int theSize = 0): capacity(max(theSize,
        int(MIN_CAPACITY))), size(theSize), items(rawMemory<ITEM>(capacity))
    {//array 0 and size > 0 means set capacity to size
        assert(size >= 0);
        if(array) for(int i = 0; i < size; ++i) new(&items[i])ITEM(array[i]);
    }
    explicit Vector(int initialSize, ITEM const& value = ITEM()): size(0),
        capacity(max(initialSize, int(MIN_CAPACITY))),
        items(rawMemory<ITEM>(capacity))
        {for(int i = 0; i < initialSize; ++i) append(value);}
    Vector(Vector const& rhs): capacity(rhs.capacity), size(rhs.size),
        items(rawMemory<ITEM>(capacity))
        {for(int i = 0; i < size; ++i) new(&items[i])ITEM(rhs.items[i]);}
    Vector& operator=(Vector const& rhs){return genericAssign(*this, rhs);}
    ~Vector(){rawDestruct(items, size);}

    void append(ITEM const& item)
    {
        if(size >= capacity) resize();
        new(&items[size++])ITEM(item);
    }
    void removeLast()
    {
        assert(size > 0);
        items[--size].~ITEM();
        if(capacity > MIN_CAPACITY && size * 4 < capacity) resize();
    }
    void swapWith(Vector& other)
    {
        swap(items, other.items);
        swap(size, other.size);
        swap(capacity, other.capacity);
    }

    ITEM const& lastItem()const{return items[size - 1];}
    ITEM& lastItem(){return items[size - 1];}
    void reverse(int left, int right)
        {while(left < right) swap(items[left++], items[right--]);}
    void reverse(){reverse(0, size - 1);}
    void appendVector(Vector const& rhs)
        {for(int i = 0; i < rhs.getSize(); ++i) append(rhs[i]);}
    Vector& operator+=(Vector const& rhs)
    {
        assert(size == rhs.size);
        for(int i = 0; i < size; ++i) items[i] += rhs.items[i];
        return *this;
    }
    Vector& operator-=(Vector const& rhs)
    {
        assert(size == rhs.size);
        for(int i = 0; i < size; ++i) items[i] -= rhs.items[i];
        return *this;
    }
    template<typename SCALAR> Vector& operator*=(SCALAR const& scalar)
    {
        for(int i = 0; i < size; ++i) items[i] *= scalar;
        return *this;
    }
    friend Vector operator+(Vector const& a, Vector const& b)
    {
        Vector result(a);
        return result += b;
    }
    friend Vector operator-(Vector const& a, Vector const& b)
    {
        Vector result(a);
        return result -= b;
    }
    friend Vector operator*(Vector const& a, ITEM const& scalar)
    {
        Vector result(a);
        return result *= scalar;
    }
    friend ITEM operator*(Vector const& a, Vector const& b)
    {
        assert(a.size == b.size);
        ITEM result(0);
        for(int i = 0; i < a.size; ++i) result += a[i] * b[i];
        return result;
    }
    Vector operator-(){return *this * -1;}
    bool operator==(Vector const& rhs)const
    {
        if(size == rhs.size)
        {
            for(int i = 0; i < size; ++i)
                if(items[i] != rhs[i]) return false;
            return true;
        }
        return false;
    }
};

}
#endif
