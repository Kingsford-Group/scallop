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
	if(argc != 3)
	{
		cout<<"usage: "<<argv[0]<< " <cuff.tmap> <ref-size>"<<endl;
		return 0;
	}

	cuffroc roc(argv[1], atoi(argv[2]));
	roc.solve();

    return 0;
}
