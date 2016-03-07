#include "disjoint_sets.h"

#include <vector>
#include <cstdio>

using namespace std;

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
