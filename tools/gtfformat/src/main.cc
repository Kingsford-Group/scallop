#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>

#include "genome.h"

using namespace std;

int main(int argc, const char **argv)
{
	if(argc != 3)
	{
		cout<<"usage: "<<argv[0]<< " <in-gtf-file> <out-gtf-file>"<<endl;
		return 0;
	}

	genome g(argv[1]);
	g.write(argv[2]);

    return 0;
}
