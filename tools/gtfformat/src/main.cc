#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>

#include "genome.h"

using namespace std;

int main(int argc, const char **argv)
{
 	if(argc != 3 && argc != 4)
	{
		cout<<"usage: "<<argv[0]<< " <in-gtf-file> <out-gtf-file> [min-coverage]"<<endl;
		return 0;
	}

	genome g(argv[1]);

	if(argc == 4)
	{
		for(int i = 0; i < g.genes.size(); i++) 
		{
			g.genes[i].filter_low_coverage_transcripts(atof(argv[3]));
		}
	}

	g.write(argv[2]);

    return 0;
}
