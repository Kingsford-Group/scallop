#include "disjoint_sets.h"

#include <cstdio>

using namespace std;

vector<int> get_representatives(disjoint_sets_t &ds, int n)
{
	vector<int> v;
	for(int i = 0; i < n; i++)
	{
		if(ds.find_set(i) != i) continue;
		v.push_back(i);
	}
	return v;
}

vector< vector<int> > get_disjoint_sets(disjoint_sets_t &ds, int n)
{
	vector< vector<int> > v;
	v.resize(n);
	for(int i = 0; i < n; i++)
	{
		int r = ds.find_set(i);
		v[r].push_back(i);
	}
	return v;
}

int test_disjoint_sets()
{
	int N = 5;
	disjoint_sets_t ds(N);
	for(int i = 0; i < N; i++) ds.make_set(i);
	ds.union_set(0, 1);
	ds.union_set(2, 3);

	for(int i = 0; i < N; i++)
	{
		printf("%d -> %d\n", i, ds.find_set(i));
	}

	return 0;
}
