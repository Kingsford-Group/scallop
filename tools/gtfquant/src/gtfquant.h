#ifndef __CUFFROC_H__
#define __CUFFROC_H__

#include "quantitem.h"
#include "genome.h"
#include <vector>
#include <map>

using namespace std;

typedef pair<double, double> PDD;
typedef pair<string, PDD> PSDD;
typedef map<string, PDD> MSDD;

class gtfquant
{
public:
	gtfquant(const string &quantfile, const string &gtffile, double m1, double m2);

public:
	genome gm;
	vector<quantitem> items;
	map<string, int> t2i;		// transcript to item
	double min_tpm;
	double min_numreads;
	MSDD msd;

public:
	int read_quant(const string &file);
	int build_indices();
	int build_msd();
	int filter();
	int compare();
};

#endif
