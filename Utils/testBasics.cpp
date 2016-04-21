#include "Stack.h"
#include "Queue.h"
#include "GCFreelist.h"
using namespace igmdk;

struct Fat
{
	enum{SIZE = 10};
	int array[SIZE];
	Fat(int last)
	{
		for(int i = 0; i < SIZE-1; ++i)
		{
			array[i] = i;
		}
		array[SIZE-1] = last;
	}
	bool operator==(Fat const& rhs)const
	{
		for(int i = 0; i < SIZE; ++i)
		{
			if(array[i] != rhs.array[i]) return false;
		}
		return true;
	}
	bool operator<(Fat const& rhs)const
	{
		for(int i = 0; i < SIZE; ++i)
		{
			if(array[i] < rhs.array[i]) return true;
			if(array[i] > rhs.array[i]) return false;
		}
		return false;
	}
	int getSize()const{return SIZE;}
	int const operator[](int i)const{return array[i];}
};

void testStackQueue()
{
    int N = 100000;
	Stack<int> a;

	for(int i=0; i < N; ++i)
	{
		a.push(i);
	}
	for(int i=0; i < N; ++i)
	{
		a.pop();
	}
	for(int i=0; i < N; ++i)
	{
		a.push(i);
	}
	for(int i=0; i < N; ++i)
	{
		a.pop();
	}
	Queue<int> b;
	for(int i=0; i < N; ++i)
	{
		b.push(i);
	}
	for(int i=0; i < N; ++i)
	{
		b.pop();
	}
	for(int i=0; i < N; ++i)
	{
		b.push(i);
	}
	for(int i=0; i < N; ++i)
	{
		b.pop();
	}
}



void DDDVector()
{
    Vector<int> Vector0to4;
    for(int i = 0; i < 5; ++i ) Vector0to4.append(i);

    cout << "breakpoint" << endl;
}

/*
void DDDDeque()
{
    Deque<char> hello;
    hello.append('h');
    hello.append('e');
    hello.append('l');
    hello.append('l');
    hello.append('o');
    cout << "breakpoint" << endl;
}*/

void DDDStack()
{
    Stack<int> Stack0to4;
    for(int i = 0; i < 5; ++i ) Stack0to4.push(i);

    cout << "breakpoint" << endl;
}

void DDDQueue()
{
    Queue<int> Queue1to4;
    for(int i = 0; i < 5; ++i ) Queue1to4.push(i);
    Queue1to4.pop();

    cout << "breakpoint" << endl;
}

void DDDList()
{
    SimpleDoublyLinkedList<int> List0to2;
    typedef SimpleDoublyLinkedList<int>::Node N;
    for(int i = 2; i >= 0; --i ) List0to2.prepend(new N(i, 0));

    cout << "breakpoint" << endl;
}

void DDDFreelist()
{
    Freelist<int> Freelist0to14R5(8);
    int* items[15];
    for(int i = 0; i < 15; ++i ) items[i] = new(Freelist0to14R5.allocate())int(i);
    for(int i = 0; i < 5; ++i ) Freelist0to14R5.remove(items[i]);

    cout << "breakpoint" << endl;
}

int main()
{
    DDDVector();
    //DDDDeque();
    DDDStack();
    DDDQueue();
    DDDList();
    DDDFreelist();
    testStackQueue();
    return 0;
}
