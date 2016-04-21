
#include "File.h"
#include "EMVector.h"
#include "../Utils/Debug.h"
#include <string>
#include "../Compression/Compression.h"
using namespace std;

using namespace igmdk;

int main(int argc, char *argv[])
{
    int start = clock();
    if(argc < 5) {DEBUG("Arguments not Provided"); return 0;}
    File in(argv[1], false), out(argv[2], true);
    bool compress;
    if(string(argv[3]) == "compress") compress = true;
    else if(string(argv[3]) == "extract") compress = false;
    else{DEBUG("Action Unknown"); return 0;}

    char method;
    enum{HUF, BWT, LZW};
    if(string(argv[4]) == "huf") method = HUF;
    else if(string(argv[4]) == "bwt") method = BWT;
    else if(string(argv[4]) == "lzw") method = LZW;
    else{DEBUG("Method Unknown"); return 0;}

    enum{N = 8096};
    char buffer[N];
    Vector<unsigned char> original, v;
    for(;;)
    {
        int size = min<long long>(N, in.bytesLeft());
        in.read(buffer, size);
        for(int i = 0; i < size; ++i)
        {
            original.append(buffer[i]);
        }
        if(size < N) break;
    }
    if(compress)
    {
        if(method == LZW)
        {
            BitStream result;
            BitStream in(original);
            LZWCompress(in, 16, result);
            v = ExtraBitsCompress(result.bitset);
        }
        else if(method == BWT)
        {
            v = BWTCompress(original);
        }
        else if(method == HUF)
        {
            v = HuffmanCompress(original);
        }
    }
    else
    {
        if(method == LZW)
        {
            BitStream in(ExtraBitsUncompress(original));
            BitStream result;
            LZWUncompress(in, 16, result);
            v = result.bitset.getStorage();
        }
        else if(method == BWT)
        {
            v = BWTUncompress(original);
        }
        else if(method == HUF)
        {
            v = HuffmanUncompress(original);
        }
    }
    out.write((char*)v.getArray(), v.getSize());
    DEBUG(clock()-start);
    return 0;
}
