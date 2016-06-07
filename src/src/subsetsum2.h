#ifndef __SUBSETSUM_H__
#define __SUBSETSUM_H__

#include "equation.h"
#include <vector>

using namespace std;

typedef pair<int, int> PI;

#define SMIN 0.0000000001

class subsetsum2
{
public:
	subsetsum2(const vector<int> &v);

private:
	const vector<int> &raw;				// raw data from input
	vector<PI> seeds;					// seeds used for sorting and rescaling
	int ubound;							// upper bound we would try 

	vector< vector<int> > table;		// dp table

public:
	vector<equation> eqns;				// multiple solutions

public:
	int solve();

private:
	int init_seeds();
	int rescale();
	int init_table();
	int fill_table();
	int optimize();
	int backtrace(int opt, vector<int> &subset);
	int recover(int &opt, vector<int> &subset);
};

#endif
