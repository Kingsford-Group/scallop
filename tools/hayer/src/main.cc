#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>

#include "hayer.h"

using namespace std;

int main(int argc, const char **argv)
{
	if(argc != 4)
	{
		cout<<"usage: "<<argv[0]<< " <input-hayer-file> <output-gtf-file> <T1|ER|EP>"<<endl;
		return 0;
	}

	hayer h;
	h.process(argv[1], argv[3]);
	h.write(argv[2]);

    return 0;
}
