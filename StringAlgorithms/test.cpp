#include <string>
#include <iostream>
#include <cstdlib>
#include <ctime>
using namespace std;
#include "StringAlgs.h"
using namespace igmdk;


int find(string const& text, string const& pattern)
{
	//Horspool matcher((unsigned char*)text.c_str(), text.length(), (unsigned char*)pattern.c_str(), pattern.length());
	//ShiftAnd matcher((unsigned char*)text.c_str(), text.length(), (unsigned char*)pattern.c_str(), pattern.length());
	HashQ<unsigned char*, Q3Hash> matcher((unsigned char*)text.c_str(), text.length(), (unsigned char*)pattern.c_str(), pattern.length());
	return matcher.findNext();
}

void timeSRT()
{
	string pattern = "hiho";
	string text = "iho 4 score and seven years ago our fathers brought forth and this continent a new hiho blah blah blah blahiho";
	int result = find(text, pattern);
	cout << result << endl;
	if(result != -1)
	{
		cout << text[result] << endl;
	}
}

void testWuManber()
{
	string pattern = "hiho";
	string pattern2 = "forth";
	string text = "iho 4 score and seven years ago our fathers brought forth and this continent a new hiho blah blah blah blahiho";
	Vector<pair<unsigned char*, int> > patterns;
    patterns.append(pair<unsigned char*, int>((unsigned char*)pattern.c_str(), pattern.length()));
    patterns.append(pair<unsigned char*, int>((unsigned char*)pattern2.c_str(), pattern2.length()));
    Vector<int> results;
    WuManber<> matcher((unsigned char*)text.c_str(), text.length(), patterns);

    int position = matcher.findNext(results);
    DEBUG(results.getSize());
    if(results.getSize() > 0)
    {
        DEBUG(position);
        for(int i = 0; i < results.getSize(); ++i)
        {
            DEBUG(results[i]);
        }
    }
}

void testLCS()
{
    Vector<unsigned char> x, y;
    x.append('h');
    x.append('e');
    x.append('h');
    x.append('u');

    y.append('h');
    y.append('i');
    y.append('h');
    Vector<Edit> result = Diff<unsigned char>(x, y);
    for(int i = 0; i < result.getSize(); ++i)
    {
        DEBUG(result[i].isInsert);
        DEBUG(result[i].position);
    }

}

void testRE()
{
    RegularExpressionMatcher re("((A*B|AC)D)");

    cout << "breakpoint" << endl;
    DEBUG(re.matches("ABCCBD"));
    DEBUG(re.matches("BCD"));
    DEBUG(re.matches("ABD"));
    DEBUG(re.matches("ACD"));
    DEBUG(re.matches("AABD"));
}



void DDDLCS()
{
    Vector<unsigned char> x, y;
    x.append('s');
    x.append('i');
    x.append('n');
    x.append('k');

    y.append('t');
    y.append('h');
    y.append('i');
    y.append('n');
    y.append('k');
    Vector<Edit> SinkIntoThink = Diff<unsigned char>(x, y);

    cout << "breakpoint" << endl;

}

int main()
{
    DDDLCS();
	clock_t start = clock();
	testRE();
	int tFL = (clock() - start);
    cout << "FL: "<<tFL << endl;
    testLCS();

	timeSRT();

	for(int i = 0; i < 1; ++i)
	{
		testWuManber();
	}
	return 0;
}
