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
		cout<<"usage: "<<argv[0]<< " <command> <gtf-file-1> <gtf-file-2>"<<endl;
		return 0;
	}

	parse_parameters(argc, argv);

	if(string(argv[1]) == "union")
	{
		gtfmerge gm;
		for(int k = 2; k < argc; k++)
		{
			gm.add_genome(argv[k]);
		}
		gm.print();
	}

    return 0;
}
