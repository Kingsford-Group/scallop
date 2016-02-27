#include "subsetsum.h"
#include <cstdio>
#include <cmath>
#include <climits>
#include <algorithm>

subsetsum::subsetsum(const vector<int> &_s, const vector<int> &_t)
	: s(_s), t(_t)
{
	augment();
	compute_closest_pair();
}

int subsetsum::augment()
{
	ss = s;
	vector<int> vi;
	int start = 0;
	for(int i = 0; i < s.size(); i++) vi.push_back(i + 1);
	for(int k = 0; k < s.size() - 1; k++) augment(s, ss, vi, start);

	tt = t;
	vi.clear();
	start = 0;
	for(int i = 0; i < t.size(); i++) vi.push_back(i + 1);
	for(int k = 0; k < t.size() - 1; k++) augment(t, tt, vi, start);

	return 0;
}

int subsetsum::augment(const vector<int> &r, vector<int> &x, vector<int> &vi, int &start)
{
	int n = x.size();
	for(int i = start; i < n; i++)
	{
		for(int j = vi[i]; j < r.size(); j++)
		{
			//printf("push back (%5d, %5d) -> (%5d, %5d)\n", i, j, s[i] + s[j], j + 1);
			x.push_back(x[i] + r[j]);
			vi.push_back(j + 1);
		}
	}
	start = n;
	return 0;
}

int subsetsum::print()
{
	printf("input first:\n");
	for(int i = 0; i < s.size(); i++) printf("%5d ", s[i]);
	printf("\n");
	printf("input second:\n");
	for(int i = 0; i < t.size(); i++) printf("%5d ", t[i]);
	printf("\n");

	printf("augment first:\n");
	for(int i = 0; i < ss.size(); i++)
	{
		if(i > 0 && i % 30 == 0) printf("\n");
		printf("%5d ", ss[i]);
	}
	printf("\n");

	printf("augment second:\n");
	for(int i = 0; i < tt.size(); i++)
	{
		if(i > 0 && i % 30 == 0) printf("\n");
		printf("%5d ", tt[i]);
	}
	printf("\n");

	printf("distance = %5d\n", dist);

	return 0;
}

int subsetsum::compute_closest_pair()
{
	sort(ss.begin(), ss.end());
	sort(tt.begin(), tt.end());
	int xi = 0;
	int yi = 0;
	dist = INT_MAX;

	while(xi < ss.size() && yi < tt.size())
	{
		int d = (int)(fabs(ss[xi] - tt[yi]));
		if(d < dist) dist = d;
		if(dist == 0) break;
		if(ss[xi] < tt[yi]) xi++;
		else yi++;
	}
	return 0;
}
