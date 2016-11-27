#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>

#include "genome.h"
#include "compare.h"
#include "config.h"

using namespace std;

int main(int argc, const char **argv)
{
	if(argc < 3)
	{
		cout<<"usage: "<<argv[0]<< " <gtf-file-1> <gtf-file-2> [-l length -a algo -s]"<<endl;
		return 0;
	}

	parse_parameters(argc, argv);

	genome g1(argv[1]);
	genome g2(argv[2]);

	if(algo == 1) compare_genome1(g1, g2);
	if(algo == 2) compare_genome2(g1, g2);

    return 0;
}
