#ifndef __CUFFROC_H__
#define __CUFFROC_H__

#include "cuffitem.h"
#include "genome.h"
#include <vector>
#include <map>

using namespace std;

class gtfcuff
{
public:
	gtfcuff(const string &cufffile, const string &gtffile);

public:
	genome gm;
	vector<cuffitem> items;
	map<string, int> t2i;		// transcript to item-index
	map<string, int> t2e;		// transcript to num-exon
	map<string, char> t2s;		// transcript to strand
	int refsize;

public:
	int read_cuff(const string &file);
	int build_indices();
	int filter_items();
	int classify(const string &f1, const string &f2);
	int roc();
	int print();
};

#endif
