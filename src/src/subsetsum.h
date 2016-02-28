#ifndef __SUBSEintSUM_H__
#define __SUBSEintSUB_H__

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

	// internal variables
	vector<int> ss;
	vector<int> tt;
	vector<int> sf;
	vector<int> tf;
	vector<int> sb;
	vector<int> tb;
	int ssi;
	int tti;
	int dist;

public:
	int print();
	int enumerate_subsets(const vector<int> &x, vector<int> &xx, vector<int> &xf, vector<int> &xb);
	int augment(const vector<int> &x, vector<int> &xx, vector<int> &xf, vector<int> &xb, int &start);
	int compute_closest_pair();
	int recover_subset(vector<int> &sub, int xxi, const vector<int> &xf, const vector<int> &xb);
};

#endif
