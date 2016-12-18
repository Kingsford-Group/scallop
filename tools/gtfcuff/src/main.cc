#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <cassert>

#include "item.h"
#include "gtfcuff.h"

using namespace std;

int main(int argc, const char **argv)
{
	if(argc != 4 && argc != 6)
	{
		cout<<"usage: " <<endl;
		cout<<"       " <<argv[0] << " roc <cuff.tmap> <ref-size>"<<endl;
		cout<<"       " <<argv[0] << " classify <cuff.tmap> <pred-file> <true-file> <false-file> "<<endl;
		return 0;
	}

	if(string(argv[1]) == "roc")
	{
		assert(argc == 4);
		gtfcuff cuff(argv[2]);
		cuff.roc(atoi(argv[3]));
	}

	if(string(argv[1]) == "classify")
	{
		assert(argc == 6);
		gtfcuff cuff(argv[2]);
		cuff.assign_pred(argv[3]);
		cuff.classify(argv[4], argv[5]);
	}

    return 0;
}
