
#include "File.h"
#include "EMVector.h"
#include "EMBTree.h"
#include "../Utils/Debug.h"
using namespace std;

using namespace igmdk;

void DDDEMVector()
{
    EMVector<int> EMVector0to9B16("EMVector.igmdk", 16);
    int K = 10;
    for(int i = 0; i < K; ++i)
    {
        EMVector0to9B16.append(i);
    }
    cout << "breakpoint" << endl;
}

int main()
{
    DDDEMVector();
    {
        EMBPlusTree<int, int> trie("124.igmdk","125.igmdk");
        int N = 15000;
        for(int i = 0; i < N; ++i)
        {
            trie.insert(-i, -i);
        }
        DEBUG("Done inserting");
        for(int j = 0; j < 1; ++j)
        {
            for(int i = 0; i < N; ++i)
            {
                bool status;
                int item = trie.find(-i, status);
                assert(status);
                assert(item == -i);
                trie.remove(-i);
            }
        }
        DEBUG(trie.getRoot());
    }
        File::remove("124.igmdk");
        File::remove("125.igmdk");
    {
        EMVector<int> vec("111.igmdk");
        int K = 1500000;
        for(int i = 0; i < K; ++i)
        {

            vec.append(-i);
        }
        DEBUG(vec.getSize());
        IOSort(vec);
        for(int i = 0; i < vec.getSize(); ++i)
        {
            //DEBUG(vec[i]);
        }
    }
    File::remove("111.igmdk");

	return 0;
}
