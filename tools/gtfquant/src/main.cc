#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <cassert>

#include "genome.h"
#include "gtfquant.h"

using namespace std;

int main(int argc, const char **argv)
{
	if(argc != 7 && argc != 6)
	{
		cout<<"usage: "<<argv[0]<< " <command> <salmon.quan> <gtffile> [min-tpm] [min-numreads] [TPM/RPKM]"<<endl;
		return 0;
	}

	gtfquant roc(argv[2], argv[3], atof(argv[4]), atof(argv[5]));

	if(string(argv[1]) == "filter") 
	{
		assert(argc == 6);
		roc.filter();
		return 0;
	}

	if(string(argv[1]) == "compare")
	{
		assert(argc == 7);
		if(string(argv[6]) == "RPKM") roc.gm.assign_TPM_by_RPKM();
		roc.compare();
	}

    return 0;
}
