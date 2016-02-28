#include "subsetsum.h"
#include <cstdio>
#include <cmath>
#include <climits>
#include <algorithm>

subsetsum::subsetsum(const vector<int> &_s, const vector<int> &_t)
	: s(_s), t(_t)
{
	enumerate_subsets(s, ss);
	enumerate_subsets(t, tt);
	compute_closest_pair();
}

int subsetsum::enumerate_subsets(const vector<int> &x, vector<int> &xx)
{
	xx = x;
	vector<int> vi;
	int start = 0;
	for(int i = 0; i < x.size(); i++) vi.push_back(i + 1);
	for(int k = 0; k < x.size() - 1; k++) augment(x, xx, vi, start);
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

	printf("distance = %5d from (%5d, %5d)\n", dist, ssi, tti);

	return 0;
}

int subsetsum::compute_closest_pair()
{
	typedef pair<int, int> PI;
	vector<PI> sss;
	vector<PI> ttt;
	for(int i = 0; i < ss.size(); i++) sss.push_back(PI(ss[i], i));
	for(int i = 0; i < tt.size(); i++) ttt.push_back(PI(tt[i], i));

	sort(sss.begin(), sss.end());
	sort(ttt.begin(), ttt.end());

	int si = 0;
	int ti = 0;
	int dist = INT_MAX;

	while(si < sss.size() && ti < ttt.size())
	{
		int d = (int)(fabs(sss[si].first - ttt[ti].first));
		if(d < dist)
		{
			dist = d;
			ssi = sss[si].second;
			tti = ttt[ti].second;
		}

		if(dist == 0) break;
		if(sss[si].first < ttt[ti].first) si++;
		else ti++;
	}

	return 0;
}
