/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "util.h"

vector<int> get_random_permutation(int n)
{
	vector<int> v;
	for(int i = 0; i < n; i++) v.push_back(i);
	for(int i = 0; i < n; i++)
	{
		int k = rand() % (n - i);
		int x = v[k];
		v[k] = v[n - i - 1];
		v[n - i - 1] = x;
	}
	return v;
}


