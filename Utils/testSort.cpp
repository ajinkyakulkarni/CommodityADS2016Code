#include "Sort.h"
#include "Vector.h"
#include <iostream>
#include <cstdlib>
using namespace igmdk;

struct Fat
{
	int key;
	int fat[100];
	bool operator<(Fat const& rhs)
	{
		return key < rhs.key;
	}
	bool operator==(Fat const& rhs)
	{
		return key == rhs.key;
	}
	Fat(int theKey):key(theKey){}
};
void timeSRT2()
{
    int simple [] = {5, 1, 3 ,0,4,2};
    Stack<int> st;
    st.push(5);
    for(int i = 0; i < 6; ++i)
    {
        DEBUG(incrementalQuickSelect((int*)simple, i, st, DefaultComparator<int>()));
    }
}

void timeSRT()
{
	Vector<int> v;
	int N = 20000;
	for(int i = N; i >= 0; --i)
	{
		v.append(rand());
	}
	bool* selected = new bool[v.getSize()];
	for(int i = 0; i < N; ++i) selected[i] = 0;
	selected[2] = true;
	selected[3] = true;
	multipleQuickSelect(v.getArray(), selected, 0, v.getSize()-1, DefaultComparator<int>());
	cout << v[2] << endl;
	cout << v[2] << endl;
    quickSort(v.getArray(), 0, v.getSize()-1, DefaultComparator<int>());
	mergeSort(v.getArray(), v.getSize(), DefaultComparator<int>());
	KSort(v.getArray(), v.getSize(), 35000, IdentityHash<int>());
}

void testMultikey()
{
    string s[7];
    s[0] = "fsdlfjl";
    s[1] = "wejk";
    s[2] = "iosufrwhrew";
    s[3] = "wqjklhdsaiohd";
    s[4] = "wioeurksd";
    s[5] = "";
    s[6] = "w";
    multikeyQuicksortNR(s, 0, 6, VectorComparator<string>());
    for(int i = 0; i < 7 ; ++i)
    {
        cout << s[i] << endl;
    }
    int k = 3;
    multikeyQuickselect(s, 0, 6, k, VectorComparator<string>());
}

int main()
{
	clock_t start = clock();
	srand(time(0));
	timeSRT();
	timeSRT2();
	testMultikey();
	int tFL = (clock() - start);
    cout << "FL: "<<tFL << endl;
	system("PAUSE");
	return 0;
}
