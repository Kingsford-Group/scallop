#ifndef __SUBSETSUM_H__
#define __SUBSETSUM_H__

#include <vector>

using namespace std;

typedef pair<int, int> PI;

class subsetsum
{
public:
	subsetsum(const vector<int> &v);

private:
	const vector<int> &raw;				// raw data from input
	vector<PI> seeds;					// seeds used for sorting and rescaling
	int target;							// target of summing up
	int ubound;							// upper bound we would try 

	vector< vector<bool> > table;		// dp table
	vector< vector<bool> > btptr;		// backtrace pointer (true means using the current number)

public:
	vector< vector<int> > subsets;		// multiple solutions
	vector<int> opts;					// errors to target

public:
	int solve();
	static int test();

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
