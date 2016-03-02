#ifndef __SUBSETSUM_H__
#define __SUBSETSUM_H__

#include <vector>

using namespace std;

class subsetsum
{
public:
	subsetsum(const vector<int> &s, const vector<int> &t);

public:
	// given two sets
	const vector<int> &s;
	const vector<int> &t;

	// output two subsets
	vector<int> subs;
	vector<int> subt;
	int dist;
	double ratio;

public:
	int solve();
	int print();
	static int test();

private:
	int enumerate_subsets(const vector<int> &x, vector<int> &xx, vector<int> &xf, vector<int> &xb);
	int augment(const vector<int> &x, vector<int> &xx, vector<int> &xf, vector<int> &xb, int &start);
	int compute_closest_pair(int &ssi, int &tti, const vector<int> &ss, const vector<int> &tt);
	int recover_subset(vector<int> &sub, int xxi, const vector<int> &xf, const vector<int> &xb);
	int compute_ratio();

};

#endif
