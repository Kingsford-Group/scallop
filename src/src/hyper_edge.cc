#include "hyper_edge.h"
#include "util.h"
#include <cstdio>
#include <algorithm>

hyper_edge::hyper_edge(const vector<int> &_v, int c)
{
	v = _v;
	count = c;
	sort(v.begin(), v.end());
}

hyper_edge::hyper_edge(const vector<int> &_v)
{
	count = 0;
	v = _v;
	sort(v.begin(), v.end());
}

bool hyper_edge::operator<(const hyper_edge &h) const
{
	if(v.size() < h.v.size()) return true;
	if(v.size() > h.v.size()) return false;

	for(int i = 0; i < v.size(); i++)
	{
		if(v[i] < h.v[i]) return true;
		if(v[i] > h.v[i]) return false;
	}
	return false;
}

int hyper_edge::increase()
{
	for(int i = 0; i < v.size(); i++) v[i]++;
	return 0;
}

int hyper_edge::print(int index) const
{
	printf("hyper edge %d: count = %d, vertices = ", index, count);
	for(int i = 0; i < v.size(); i++) printf("%d ", v[i]);
	printf("\n");
	return 0;
}
