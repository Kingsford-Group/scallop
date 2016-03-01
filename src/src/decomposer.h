#ifndef __DECOMPOSER_H__
#define __DECOMPOSER_H__

#include <vector>

using namespace std;

typedef pair<int, int> PI;
typedef pair<PI, int> PPII;
typedef vector<PPII> VPPII;
typedef pair< vector<int>, vector<int> > PVV;

class decomposer 
{
public:
	decomposer(const vector<int> &s, const vector<int> &t);

public:
	// given two sets
	const vector<int> &s;
	const vector<int> &t;

	// decomposition
	vector<PVV> subsets;
	vector<VPPII> vpi;

public:
	int solve();
	int print();
	static int test();

private:
	int build_subsets();
	vector<int> complement(const vector<int> &x, int n);
	vector<int> decipher(const vector<int> &x, const vector<int> &r);
	vector<VPPII> enumerate_paths(const vector<int> &x, const vector<int> &y);
};

#endif
