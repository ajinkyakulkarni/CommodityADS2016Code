#ifndef SORT_H
#define SORT_H
#include "Utils.h"
#include "Stack.h"
#include "../RandomNumberGeneration/Random.h"
namespace igmdk{

template<typename ITEM, typename COMPARATOR>
void insertionSort(ITEM* vector, int left, int right, COMPARATOR const& c)
{
    for(int i = left + 1; i <= right; ++i)
    {
        ITEM e = vector[i];
        int j = i;
        for(;j > left && c.isLess(e, vector[j - 1]); --j)
            vector[j] = vector[j - 1];
        vector[j] = e;
    }
}

template<typename ITEM, typename COMPARATOR>
int pickPivot(ITEM* vector, int left, int right, COMPARATOR const& c)
{
    int i = GlobalRNG.inRange(left, right), j =
        GlobalRNG.inRange(left, right), k = GlobalRNG.inRange(left, right);
    if(c.isLess(vector[j], vector[i])) swap(i, j);
    //i <= j, decide where k goes
    return c.isLess(vector[k], vector[i]) ?
        i : c.isLess(vector[k], vector[j]) ? k : j;
}

template<typename ITEM, typename COMPARATOR> void partition3(ITEM* vector,
    int left, int right, int& i, int& j, COMPARATOR const& c)
{
    ITEM p = vector[pickPivot(vector, left, right, c)];
    int lastLeftEqual = i = left - 1, firstRightEqual = j = right + 1;
    for(;;)//the pivot is the sentinel for the first pass
    {//after one swap swapped items act as sentinels
        while(c.isLess(vector[++i], p));
        while(c.isLess(p, vector[--j]));
        if(i >= j) break;
        swap(vector[i], vector[j]);
        //swap equal items to the sides
        if(c.isEqual(vector[i], p)) swap(vector[++lastLeftEqual], vector[i]);
        if(c.isEqual(vector[j], p))
            swap(vector[--firstRightEqual], vector[j]);
    }
    //invariant: i == j if they stop at an item = pivot
    //and this can happen at both left and right item
    //or they cross over and i = j + 1
    if(i == j){++i; --j;}
    //swap side items to the middle
    for(int k = left; k <= lastLeftEqual; ++k) swap(vector[k], vector[j--]);
    for(int k = right; k >= firstRightEqual; --k)
        swap(vector[k], vector[i++]);
}

template<typename ITEM, typename COMPARATOR> ITEM quickSelect(ITEM* vector,
    int left, int right, int k, COMPARATOR const& c)
{
    assert(k >= left && k <= right);
    for(int i, j; left < right;)
    {
        partition3(vector, left, right, i, j, c);
        if(k >= i) left = i;
        else if(k <= j) right = j;
        else break;
    }
    return vector[k];
}
template<typename ITEM> ITEM quickSelect(ITEM* vector, int size, int k)
    {return quickSelect(vector, 0, size - 1, k, DefaultComparator<ITEM>());}

template<typename ITEM, typename COMPARATOR> ITEM incrementalQuickSelect(
    ITEM* vector, int left, Stack<int>& s, COMPARATOR const& c)
{
    for(int right, i, j; left < (right = s.getTop()); s.push(j))
        partition3(vector, left, right, i, j, c);
    s.pop();
    return vector[left];
}
template<typename ITEM, typename COMPARATOR> void incrementalSort(
    ITEM* vector, int n, COMPARATOR const& c)
{
    Stack<int> s;
    s.push(n - 1);
    for(int i = 0; i < n; ++i) incrementalQuickSelect(vector, i, s, c);
}

template<typename ITEM, typename COMPARATOR> void multipleQuickSelect(ITEM*
    vector, bool* selected, int left, int right, COMPARATOR const& c)
{
    while(right - left > 16)
    {
        int i, j;
        for(i = left; i <= right && !selected[i]; ++i);
        if(i == right+1) return;//none are selected
        partition3(vector, left, right, i, j, c);
        if(j - left < right - i)//smaller first
        {
            multipleQuickSelect(vector, selected, left, j, c);
            left = i;
        }
        else
        {
            multipleQuickSelect(vector, selected, i, right, c);
            right = j;
        }
    }
    insertionSort(vector, left, right, c);
}

template<typename ITEM> void quickSort(ITEM* vector, int left, int right)
    {quickSort(vector, left, right, DefaultComparator<ITEM>());}
template<typename ITEM, typename COMPARATOR>
void quickSort(ITEM* vector, int left, int right, COMPARATOR const& c)
{
    while(right - left > 16)
    {
        int i, j;
        partition3(vector, left, right, i, j, c);
        if(j - left < right - i)//smaller first
        {
            quickSort(vector, left, j, c);
            left = i;
        }
        else
        {
            quickSort(vector, i, right, c);
            right = j;
        }
    }
    insertionSort(vector, left, right, c);
}

template<typename ITEM, typename COMPARATOR> void merge(ITEM* vector,
    int left, int middle, int right, COMPARATOR const& c, ITEM* storage)
{
    for(int i = left, j = middle + 1; left <= right; ++left)
    {
        bool useRight = i > middle || (j <= right &&
            c.isLess(storage[j], storage[i]));
        vector[left] = storage[(useRight ? j : i)++];
    }
}
template<typename ITEM, typename COMPARATOR> void mergeSortHelper(
    ITEM* vector, int left, int right, COMPARATOR const& c, ITEM* storage)
{
    if(right - left > 16)
    {//sort storage using vector as storage
        int middle = (right + left)/2;
        mergeSortHelper(storage, left, middle, c, vector);
        mergeSortHelper(storage, middle + 1, right, c, vector);
        merge(vector, left, middle, right, c, storage);
    }
    else insertionSort(vector, left, right, c);
}
template<typename ITEM, typename COMPARATOR>
void mergeSort(ITEM* vector, int n, COMPARATOR const& c)
{
    if(n <= 1) return;
    Vector<ITEM> storage(vector, n);
    mergeSortHelper(vector, 0, n - 1, c, storage.getArray());
}
template<typename ITEM> void mergeSort(ITEM* vector, int n)
    {mergeSort(vector, n, DefaultComparator<ITEM>());}

template<typename VECTOR, typename COMPARATOR> void multikeyQuicksort(VECTOR*
    vector, int left, int right, COMPARATOR const& c)
{
    if(right - left < 1) return;
    int i, j;
    partition3(vector, left, right, i, j, c);
    ++c.depth;
    multikeyQuicksort(vector, j + 1, i - 1, c);
    --c.depth;
    multikeyQuicksort(vector, left, j, c);
    multikeyQuicksort(vector, i, right, c);
}

template<typename VECTOR, typename COMPARATOR> void multikeyQuicksortNR(
    VECTOR* vector, int left, int right, COMPARATOR const& c,
    int maxDepth = numeric_limits<int>::max())
{
    Stack<int> stack;
    stack.push(left);
    stack.push(right);
    stack.push(0);
    while(!stack.isEmpty())
    {
        c.depth = stack.pop();
        right = stack.pop();
        left = stack.pop();
        if(right - left > 0 && c.depth < maxDepth)
        {
            int i, j;
            partition3(vector, left, right, i, j, c);
            //left
            stack.push(left);
            stack.push(j);
            stack.push(c.depth);
            //right
            stack.push(i);
            stack.push(right);
            stack.push(c.depth);
            //middle
            stack.push(j + 1);
            stack.push(i - 1);
            stack.push(c.depth + 1);
        }
    }
}

template<typename VECTOR, typename COMPARATOR> void multikeyQuickselect(
    VECTOR* vector, int left, int right, int k, COMPARATOR const& c)
{
    assert(k >= left && k <= right);
    for(int d = 0, i, j; right - left >= 1;)
    {
        partition3(vector, left, right, i, j, c);
        if(k <= j) right = j;
        else if (k < i)
        {
            left = j + 1;
            right = i - 1;
            ++c.depth;
        }
        else left = i;
    }
}

template<typename VECTOR> struct VectorComparator
{
    mutable int depth;
    VectorComparator(): depth(0){}
    bool isLess(VECTOR const& lhs, VECTOR const& rhs)const
    {
        return depth < lhs.length() ?
            depth < rhs.length() && lhs[depth] < rhs[depth] :
            depth < rhs.length();
    }
    bool isEqual(VECTOR const& lhs, VECTOR const& rhs)const
    {
        return depth < lhs.length() ?
            depth < rhs.length() && lhs[depth] == rhs[depth] :
            depth >= rhs.length();
    }
};

template<typename ITEM, typename COMPARATOR> int binarySearch(ITEM const*
    vector, int left, int right, ITEM const& key, COMPARATOR const& c)
{
    while(left <= right)
    {
        int middle = (left + right)/2;
        if(c.isEqual(key, vector[middle])) return middle;
        c.isLess(key, vector[middle]) ?
            right = middle - 1 : left = middle + 1;
    }
    return -1;
}

template<typename ITEM> void permutationSort(ITEM* a, int* permutation,
    int n)
{
    for(int i = 0; i < n; ++i)
        if(permutation[i] != i)
        {
            ITEM temp = a[i];
            int from = i, to;
            for(;;)
            {
                from = permutation[to = from];
                permutation[to] = to;//mark processed
                if(from == i) break;
                a[to] = a[from];
            }
            a[to] = temp;//complete cycle
        }
}

void countingSort(int* vector, int n, int N)
{
    Vector<int> counter(N, 0);
    for(int i = 0; i < n; ++i) ++counter[vector[i]];
    for(int i = 0, index = 0; i < N; ++i)
        while(counter[i]-- > 0) vector[index++] = i;
}

template<typename ITEM> struct IdentityHash
    {int operator()(ITEM const& item)const{return item;}};

template<typename ITEM, typename ORDERED_HASH> void KSort(ITEM* a, int n,
    int N, ORDERED_HASH const& h)
{
    ITEM* temp = rawMemory<ITEM>(n);
    Vector<int> count(N + 1, 0);
    for(int i = 0; i < n; ++i) ++count[h(a[i]) + 1];
    for(int i = 0; i < N; ++i) count[i + 1] += count[i];//accumulate counts
    //rearrange items
    for(int i = 0; i < n; ++i) new(&temp[count[h(a[i])]++])ITEM(a[i]);
    for(int i = 0; i < n; ++i) a[i] = temp[i];
    rawDelete(temp);
}
}//end namespace
#endif
