#ifndef COMPRESSION_H
#define COMPRESSION_H
#include "../RandomTreap/Trie.h"
#include "../Heaps/Heap.h"
#include "../StringAlgorithms/SuffixArray.h"
#include "Stream.h"
#include <cstdlib>
namespace igmdk{

Vector<unsigned char> ExtraBitsCompress(Bitset<unsigned char> bitset)
{
    bitset.getStorage().append(bitset.garbageBits());
    return bitset.getStorage();
}
Bitset<unsigned char> ExtraBitsUncompress(Vector<unsigned char> byteArray)
{
    int garbageBits = byteArray.lastItem();
    byteArray.removeLast();
    Bitset<unsigned char> result(byteArray);
    while(garbageBits--) result.removeLast();
    return result;
}

void byteEncode(unsigned long long n, BitStream& result)
{
    do
    {
        unsigned char r = n % 128;
        n /= 128;
        if(n) r += 128;
        result.writeByte(r);
    }while(n);
}
unsigned long long byteDecode(BitStream& stream)
{
    unsigned long long n = 0, base = 1;
    for(;; base *= 128)
    {
        unsigned char code = stream.readByte(), value = code % 128;
        n += base * value;
        if(value == code) break;
    }
    return n;
}

void UnaryEncode(int n, BitStream& result)
{
    while(n--) result.writeBit(true);
    result.writeBit(false);
}
int UnaryDecode(BitStream& code)
{
    int n = 0;
    while(code.readBit()) ++n;
    return n;
}

void GammaEncode(unsigned long long n, BitStream& result)
{
    assert(n > 0);
    int N = lgFloor(n);
    UnaryEncode(N, result);
    if(N > 0) result.writeValue(n - twoPower(N), N);
}
unsigned long long GammaDecode(BitStream& code)
{
    int N = UnaryDecode(code);
    return twoPower(N) + (N > 0 ? code.readValue(N) : 0);
}

void advanceFib(unsigned long long& f1, unsigned long long& f2)
{
    unsigned long long temp = f2;
    f2 += f1;
    f1 = temp;
}
void FibonacciEncode(unsigned long long n, BitStream& result)
{
    assert(n > 0);
    //find largest fib number f1 <= n
    unsigned long long f1 = 1, f2 = 2;
    while(f2 <= n) advanceFib(f1, f2);
    //mark the numbers from highest to lowest
    Bitset<unsigned char> reverse;
    while(f2 > 1)
    {
        reverse.append(n >= f1);
        if(n >= f1) n -= f1;
        unsigned long long temp = f1;
        f1 = f2 - f1;
        f2 = temp;
    }//change order to lowest to highest and add terminator
    reverse.reverse();
    result.bitset.appendBitset(reverse);
    result.writeBit(true);
}
unsigned long long FibonacciDecode(BitStream& code)
{
    unsigned long long n = 0, f1 = 1, f2 = 2;
    for(bool prevBit = false;; advanceFib(f1, f2))
    {//add on the next Fibonacci number until see 11
        bool bit = code.readBit();
        if(bit)
        {
            if(prevBit) break;
            n += f1;
        }
        prevBit = bit;
    }
    return n;
}

Vector<unsigned char> MoveToFrontTransform(bool compress,
    Vector<unsigned char>const& byteArray)
{
    unsigned char list[256], j, letter;
    for(int i = 0; i < sizeof(list); ++i) list[i] = i;
    Vector<unsigned char> resultArray;
    for(int i = 0; i < byteArray.getSize(); ++i)
    {
        if(compress)
        {
            j = 0;
            letter = byteArray[i];
            while(list[j] != letter) ++j;
            resultArray.append(j);
        }
        else
        {
            j = byteArray[i];
            letter = list[j];
            resultArray.append(letter);
        }
        for(; j > 0; --j) list[j] = list[j - 1];
        list[0] = letter;
    }
    return resultArray;
}

enum {RLE_E1 = 255, RLE_E2 = 254};
Vector<unsigned char> RLECompress(Vector<unsigned char>const& byteArray)
{
    Vector<unsigned char> result;
    for(int i = 0; i < byteArray.getSize();)
    {
        unsigned char letter = byteArray[i++];
        result.append(letter);
        int count = 0;
        while(count < RLE_E2 - 1 && i + count < byteArray.getSize() &&
            byteArray[i + count] == letter) ++count;
        if(count > 1 || (letter == RLE_E1 && count == 1))
        {
            result.append(RLE_E1);
            result.append(count);
            i += count;
        }
        else if(letter == RLE_E1) result.append(RLE_E2);
    }
    return result;
}
Vector<unsigned char> RLEUncompress(Vector<unsigned char>const& byteArray)
{
    Vector<unsigned char> result;
    for(int i = 0; i < byteArray.getSize();)
    {
        unsigned char letter = byteArray[i++];
        if(letter == RLE_E1 && byteArray[i] != RLE_E1)
        {
            unsigned char count = byteArray[i++];
            if(count == RLE_E2) count = 1;
            else letter = result.lastItem();//need temp if vector reallocates
            while(count--) result.append(letter);
        }
        else result.append(letter);
    }
    return result;
}

void LZWCompress(BitStream& in, int maxBits, BitStream& out)
{
    if(!in.bytesLeft()) return;
    Vector<unsigned char> word(1, 1);
    TernaryTreapTrie<int> dictionary;
    TernaryTreapTrie<int>::Handle h;
    int n = 0;
    while(n < 256)
    {
        word[0] = n;
        dictionary.insert(word.getArray(), 1, n++);
    }
    word = Vector<unsigned char>();
    do
    {
        unsigned char c = in.readByte();
        word.append(c);
        if(!dictionary.findIncremental(word.getArray(), word.getSize(), h))
        {
            out.writeValue(*dictionary.find(word.getArray(),
                word.getSize() - 1), lgCeiling(n));
            if(n < twoPower(maxBits))
                dictionary.insert(word.getArray(), word.getSize(), n++);
            word = Vector<unsigned char>();
            word.append(c);
        }
    }while(in.bytesLeft());
    out.writeValue(*dictionary.find(word.getArray(), word.getSize()),
        lgCeiling(n));
}

void LZWUncompress(BitStream& in, int maxBits, BitStream& out)
{
    int size = twoPower(maxBits), n = 0, lastIndex = -1;
    Vector<Vector<unsigned char> > dictionary(size);
    for(; n < 256; ++n) dictionary[n].append(n);
    while(in.bitsLeft())
    {
        int index = in.readValue(lastIndex == -1 ? 8 :
            min(maxBits, lgCeiling(n + 1)));
        if(lastIndex != -1 && n < size)
        {
            Vector<unsigned char> word = dictionary[lastIndex];
            word.append((index == n ? word : dictionary[index])[0]);
            dictionary[n++] = word;
        }
        for(int i = 0; i < dictionary[index].getSize(); ++i)
            out.writeByte(dictionary[index][i]);
        lastIndex = index;
    }
}

struct HuffmanTree
{
    enum{W = numeric_limits<unsigned char>::digits, N = 1 << W};
    struct Node
    {
        unsigned char letter;
        int count;
        Node *left, *right;
        Node(int theCount, Node* theLeft, Node* theRight,
             unsigned char theLetter): left(theLeft), right(theRight),
             count(theCount), letter(theLetter) {}
        bool operator<(Node const& rhs)const{return count < rhs.count;}
        void traverse(Bitset<unsigned char>* codebook,
            Bitset<unsigned char>& currentCode)
        {
            if(left)
            {
                currentCode.append(false);
                left->traverse(codebook, currentCode);
                currentCode.removeLast();
                currentCode.append(true);
                right->traverse(codebook, currentCode);
                currentCode.removeLast();
            }
            else codebook[letter] = currentCode;
        }
        void append(Bitset<unsigned char>& result)
        {
            result.append(!left);
            if(left)
            {
                left->append(result);
                right->append(result);
            }
            else result.appendValue(letter, W);
        }
    }* root;
    Freelist<Node> freelist;

    HuffmanTree(Vector<unsigned char> const& byteArray)
    {//calculate frequencies
        int counts[N];
        for(int i = 0; i < N; ++i) counts[i] = 0;
        for(int i = 0; i < byteArray.getSize(); ++i) ++counts[byteArray[i]];
        //create leaf nodes
        Heap<Node*, PointerComparator<Node*> > queue;
        for(int i = 0; i < N; ++i)
            if(counts[i] > 0) queue.insert(new(freelist.allocate())
                Node(counts[i], 0, 0, i));
        //merge leaf nodes to create the tree
        while(queue.getSize() > 1)
        {
            Node *first = queue.deleteMin(), *second = queue.deleteMin();
            queue.insert(new(freelist.allocate())
                Node(first->count + second->count, first, second, 0));
        }
        root = queue.getMin();
    }

    Node* readHuffmanTree(BitStream& text)
    {
        Node *left = 0, *right = 0;
        unsigned char letter;
        if(text.readBit()) letter = text.readValue(W);
        else
        {
            left = readHuffmanTree(text);
            right = readHuffmanTree(text);
        }
        return new(freelist.allocate())Node(0, left, right, letter);
    }
    HuffmanTree(BitStream& text){root = readHuffmanTree(text);}

    void writeTree(Bitset<unsigned char>& result){root->append(result);}
    void populateCodebook(Bitset<unsigned char>* codebook)
    {
        Bitset<unsigned char> temp;
        root->traverse(codebook, temp);
    }

    void decode(BitStream& text, Vector<unsigned char>& result)
    {
        for(Node* current = root;;
            current = text.readBit() ? current->right : current->left)
        {
            if(!current->left)
            {
                result.append(current->letter);
                current = root;
            }
            if(!text.bitsLeft()) break;
        }
    }
};

Vector<unsigned char> HuffmanCompress(Vector<unsigned char>const& byteArray)
{
    HuffmanTree tree(byteArray);
    Bitset<unsigned char> codebook[HuffmanTree::N], result;
    tree.populateCodebook(codebook);
    tree.writeTree(result);
    for(int i = 0; i < byteArray.getSize(); ++i)
        result.appendBitset(codebook[byteArray[i]]);
    return ExtraBitsCompress(result);
}

Vector<unsigned char> HuffmanUncompress(
    Vector<unsigned char>const& byteArray)
{
    BitStream text(ExtraBitsUncompress(byteArray));
    HuffmanTree tree(text);
    Vector<unsigned char> result;
    tree.decode(text, result);
    return result;
}

Vector<unsigned char> BurrowsWheelerTransform(
    Vector<unsigned char> const& byteArray)
{
    int original = 0, size = byteArray.getSize();
    Vector<int> BTWArray = suffixArray<BWTRank>(byteArray.getArray(), size);
    Vector<unsigned char> result;
    for(int i = 0; i < size; ++i)
    {
        int index = BTWArray[i];
        if(index == 0)
        {
            original = i;
            index = size;
        }
        result.append(byteArray[index - 1]);
    }//make portability assumption that 4 bytes is enough
    Vector<unsigned char> code = ReinterpretEncode(original, 4);
    for(int i = 0; i < code.getSize(); ++i) result.append(code[i]);
    return result;
}

Vector<unsigned char> BurrowsWheelerReverseTransform(
     Vector<unsigned char> const& byteArray)
{
    int counts[256], firstPositions[256],
        textSize = byteArray.getSize() - 4;
    for(int i = 0; i < 256; ++i) counts[i] = 0;
    Vector<int> ranks(textSize);
    for(int i = 0; i < textSize; ++i) ranks[i] = counts[byteArray[i]]++;
    firstPositions[0] = 0;
    for(int i = 0; i < 255; ++i)
        firstPositions[i + 1] = firstPositions[i] + counts[i];
    Vector<unsigned char> index, result(textSize);
    for(int i = 0; i < 4; ++i)
        index.append(byteArray[i + textSize]);
    for(int i = textSize - 1, ix = ReinterpretDecode(index); i >= 0; --i)
        ix = ranks[ix] + firstPositions[result[i] = byteArray[ix]];
    return result;
}

Vector<unsigned char> BWTCompress(Vector<unsigned char>const& byteArray)
{
    return HuffmanCompress(RLECompress(MoveToFrontTransform(true,
        BurrowsWheelerTransform(byteArray))));
}
Vector<unsigned char> BWTUncompress(Vector<unsigned char>const& byteArray)
{
    return BurrowsWheelerReverseTransform(MoveToFrontTransform(false,
       RLEUncompress(HuffmanUncompress(byteArray))));
}

}
#endif
