#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <cassert>

#include "gtfmerge.h"
#include "config.h"

using namespace std;

int main(int argc, const char **argv)
{
	if(argc == 1)
	{
		cout<<"usage: "<<argv[0]<< " union <input-gtf-file-list> <output-merged-gtf> [options]"<<endl;
		cout<<"usage: "<<argv[0]<< " intersection <input-gtf-file-list> <output-merged-gtf> [options]"<<endl;
		return 0;
	}

	parse_parameters(argc, argv);

	if(string(argv[1]) == "union")
	{
		assert(argc >= 4);
		gtfmerge gm;
		gm.add_genomes(argv[2]);

		genome1 g1;
		gm.build_union(g1);
		g1.write(argv[3]);
	}

	if(string(argv[1]) == "intersection")
	{
		assert(argc >= 4);
		gtfmerge gm;
		gm.add_genomes(argv[2]);

		genome1 g1;
		gm.build_pairwise_intersection(g1);
		g1.write(argv[3]);
	}

    return 0;
}
