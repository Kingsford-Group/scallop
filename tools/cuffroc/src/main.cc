#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>

#include "item.h"
#include "cuffroc.h"

using namespace std;

int main(int argc, const char **argv)
{
	if(argc != 4 && argc != 5)
	{
		cout<<"usage: "<<argv[0]<< " <cuff.tmap> <gtffile> <ref-size> [min_exon_num]"<<endl;
		return 0;
	}

	if(argc == 4)
	{
		cuffroc roc(argv[1], argv[2], atoi(argv[3]));
		roc.solve();
		roc.print();
	}

	if(argc == 5)
	{
		cuffroc roc(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));
		roc.solve();
		roc.print();
	}

    return 0;
}
