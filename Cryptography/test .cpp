#include "Cryptography.h"
#include "../Utils/Debug.h"
using namespace igmdk;
#include <iostream>
#include <ctime>
using namespace std;

void timeRT()
{
    unsigned char block[128] = "top secret info", key[128] = "123456";
    Vector<unsigned char> data, password;
    for(int i = 0; i < 128; ++i){
     data.append(block[i]);
     password.append(key[i]);
    }
    Vector<unsigned char> zap = simpleDecrypt(simpleEncrypt(data, password), password);

    cout << "breakpoint" << endl;
    for(int i = 0; i < 128; ++i) block[i] = data[i];
    DEBUG(block);
}
int main()
{
	clock_t start = clock();


	timeRT();
	int tFL = (clock() - start);
    cout << "FL: "<<tFL << endl;
	return 0;
}
