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
		cout<<"usage: "<<argv[0]<< " union <merged-gtf> <gtf-file-1> <gtf-file-2> [<gtf-file-3> ...]"<<endl;
		return 0;
	}

	parse_parameters(argc, argv);

	if(string(argv[1]) == "union")
	{
		assert(argc >= 5);
		gtfmerge gm;
		for(int k = 3; k < argc; k++)
		{
			gm.add_genome(argv[k]);
		}
		gm.print();

		genome1 g1;
		gm.build_union(g1);
		g1.print(9);
		g1.write(argv[2]);
	}

    return 0;
}
