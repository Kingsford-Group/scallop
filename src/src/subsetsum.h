#ifndef __SUBSETSUM_H__
#define __SUBSETSUM_H__

#include <vector>

using namespace std;

namespace subsetsum
{
	double compute_closest_subsets(const vector<int> &s, const vector<int> &t, vector<int> &subs, vector<int> &subt);
	int enumerate_subsets(const vector<int> &x, vector<int> &xx);
	int enumerate_subsets(const vector<int> &x, vector<int> &xx, vector<int> &xf, vector<int> &xb);
	int augment(const vector<int> &x, vector<int> &xx, vector<int> &xf, vector<int> &xb, int &start);
	int compute_closest_pair(int &ssi, int &tti, const vector<int> &ss, const vector<int> &tt);
	int recover_subset(vector<int> &sub, int xxi, const vector<int> &xf, const vector<int> &xb);
	double compute_subset_ratio(const vector<int> &s, const vector<int> &t, const vector<int> &subs, const vector<int> &subt);
}

#endif
