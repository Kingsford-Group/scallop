#ifndef __CUFFROC_H__
#define __CUFFROC_H__

#include "cuffitem.h"
#include <vector>
#include <map>

using namespace std;

class cuffroc
{
public:
	cuffroc(const string &cufffile, const string &gtffile, int r, int m = -1);

public:
	vector<cuffitem> items;
	map<string, int> t2e;
	int refsize;
	int mexons;

public:
	int read_cuff(const string &file);
	int read_gtf(const string &file);
	int filter_items();
	int solve();
	int print();
};

#endif
