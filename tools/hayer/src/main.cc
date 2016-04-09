#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>

#include "hayer.h"

using namespace std;

int main(int argc, const char **argv)
{
	if(argc != 3)
	{
		cout<<"usage: "<<argv[0]<< "<input-hayer-file> <output-gtf-file>"<<endl;
		return 0;
	}

	hayer h(argv[1]);
	h.write(argv[2]);

    return 0;
}
