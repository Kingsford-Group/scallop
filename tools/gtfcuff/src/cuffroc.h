#ifndef __CUFFROC_H__
#define __CUFFROC_H__

#include "cuffitem.h"
#include "genome.h"
#include <vector>
#include <map>

using namespace std;

class cuffroc
{
public:
	cuffroc(const string &cufffile, const string &gtffile, int r, int m = -1, int f = 1, double p = 0.2);

public:
	genome gm;
	vector<cuffitem> items;
	map<string, int> t2i;		// transcript to item-index
	map<string, int> t2e;		// transcript to num-exon
	map<string, char> t2s;		// transcript to strand
	int refsize;
	int mexons;
	int ftype;
	double pratio;

public:
	int read_cuff(const string &file);
	int build_indices();
	int filter_items();
	int solve();
	int classify();
	int print();
};

#endif
