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
	vector<int> subset;					// optimal subset
	int opt;							// error to optimal solution

public:
	int solve();
	static int test();

private:
	int init_seeds();
	int rescale();
	int init_table();
	int fill_table();
	int optimize();
	int backtrace();
	int recover();
	int print();

public:
	// static functions using enumeration 
	static double compute_closest_subsets(const vector<int> &s, const vector<int> &t, vector<int> &subs, vector<int> &subt);
	static int enumerate_subsets(const vector<int> &x, vector<int> &xx);
	static int enumerate_subsets(const vector<int> &x, vector<int> &xx, vector<int> &xf, vector<int> &xb);
	static int augment(const vector<int> &x, vector<int> &xx, vector<int> &xf, vector<int> &xb, int &start);
	static int compute_closest_pair(int &ssi, int &tti, const vector<int> &ss, const vector<int> &tt);
	static int recover_subset(vector<int> &sub, int xxi, const vector<int> &xf, const vector<int> &xb);
	static double compute_subset_ratio(const vector<int> &s, const vector<int> &t, const vector<int> &subs, const vector<int> &subt);
};

#endif
