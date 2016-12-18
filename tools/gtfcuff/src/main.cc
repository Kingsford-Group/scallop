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
	if(argc != 5)
	{
		cout<<"usage: " <<endl;
		cout<<"       " <<argv[0] << " cuff <cuff.tmap> <gtffile> <ref-size>"<<endl;
		cout<<"       " <<argv[0] << " classify <cuff.tmap> <gtffile> <true-file> <false-file> "<<endl;
		return 0;
	}

	if(string(argv[1]) == "cuff")
	{
		assert(argc == 5);
		gtfcuff cuff(argv[2], argv[3]);
		cuff.refsize = atoi(argv[4]);
		cuff.roc();
	}

	if(string(argv[1]) == "classify")
	{
		assert(argc == 6);
		gtfcuff cuff(argv[2], argv[3]);
		cuff.classify(argv[4], argv[5]);
	}

    return 0;
}
