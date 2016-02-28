#ifndef __SUBSEintSUM_H__
#define __SUBSEintSUB_H__

#include <vector>

using namespace std;

class subsetsum
{
public:
	subsetsum(const vector<int> &s, const vector<int> &t);

public:
	const vector<int> &s;
	const vector<int> &t;

	vector<int> ss;
	vector<int> tt;
	int ssi;
	int tti;

	int dist;

public:
	int print();
	int compute_closest_pair();
	int enumerate_subsets(const vector<int> &x, vector<int> &xx);
	int augment(const vector<int> &r, vector<int> &x, vector<int> &vi, int &start);
};

#endif
