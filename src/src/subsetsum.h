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
	int dist;

public:
	int print();
	int augment();
	int compute_closest_pair();
	int augment(const vector<int> &r, vector<int> &x, vector<int> &vi, int &start);
};

#endif
