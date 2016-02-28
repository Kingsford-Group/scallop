#include "subsetsum.h"
#include <cstdio>
#include <cmath>
#include <climits>
#include <algorithm>
#include <cassert>

subsetsum::subsetsum(const vector<int> &_s, const vector<int> &_t)
	: s(_s), t(_t)
{
	enumerate_subsets(s, ss, sf, sb);
	enumerate_subsets(t, tt, tf, tb);
	compute_closest_pair();
	recover_subset(subs, ssi, sf, sb);
	recover_subset(subt, tti, tf, tb);
}

int subsetsum::enumerate_subsets(const vector<int> &x, vector<int> &xx, vector<int> &xf, vector<int> &xb)
{
	xx = x;
	int start = 0;
	xf.clear();
	xb.clear();
	for(int i = 0; i < x.size(); i++) xf.push_back(i + 1);
	for(int i = 0; i < x.size(); i++) xb.push_back(-1);
	for(int k = 0; k < x.size() - 1; k++) augment(x, xx, xf, xb, start);
	return 0;
}

int subsetsum::augment(const vector<int> &x, vector<int> &xx, vector<int> &xf, vector<int> &xb, int &start)
{
	int n = xx.size();
	for(int i = start; i < n; i++)
	{
		for(int j = xf[i]; j < x.size(); j++)
		{
			xx.push_back(xx[i] + x[j]);
			xf.push_back(j + 1);
			xb.push_back(i);
		}
	}
	start = n;
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
	dist = INT_MAX;

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

int subsetsum::recover_subset(vector<int> &sub, int xxi, const vector<int> &xf, const vector<int> &xb)
{
	sub.clear();
	int x = xxi;
	while(true)
	{
		assert(x >= 1 && x <= xf.size());
		sub.push_back(xf[x] - 1);
		if(xb[x] == -1) break;
		x = xb[x];
	}
	return 0;
}

int subsetsum::print()
{
	printf("input  first: ");
	for(int i = 0; i < s.size(); i++) printf("%5d ", s[i]);
	printf("\n");
	printf("input second: ");
	for(int i = 0; i < t.size(); i++) printf("%5d ", t[i]);
	printf("\n");

	printf("augment  first: ");
	for(int i = 0; i < ss.size(); i++)
	{
		if(i > 0 && i % 20 == 0) printf("\n");
		printf("%5d ", ss[i]);
	}
	printf("\n");
	printf("forward  first: ");
	for(int i = 0; i < sf.size(); i++)
	{
		if(i > 0 && i % 20 == 0) printf("\n");
		printf("%5d ", sf[i]);
	}
	printf("\n");
	printf("backward first: ");
	for(int i = 0; i < sb.size(); i++)
	{
		if(i > 0 && i % 20 == 0) printf("\n");
		printf("%5d ", sb[i]);
	}
	printf("\n");


	printf("augment  second: ");
	for(int i = 0; i < tt.size(); i++)
	{
		if(i > 0 && i % 20 == 0) printf("\n");
		printf("%5d ", tt[i]);
	}
	printf("\n");
	printf("forward  second: ");
	for(int i = 0; i < tf.size(); i++)
	{
		if(i > 0 && i % 20 == 0) printf("\n");
		printf("%5d ", tf[i]);
	}
	printf("\n");
	printf("backward second: ");
	for(int i = 0; i < tb.size(); i++)
	{
		if(i > 0 && i % 20 == 0) printf("\n");
		printf("%5d ", tb[i]);
	}
	printf("\n");


	printf("distance = %5d from (%5d, %5d)\n", dist, ssi, tti);
	printf("first  subset: ");
	for(int i = 0; i < subs.size(); i++) printf("%3d:%5d ", subs[i], s[subs[i]]);
	printf("\n");
	printf("second subset: ");
	for(int i = 0; i < subt.size(); i++) printf("%3d:%5d ", subt[i], t[subt[i]]);
	printf("\n");

	return 0;
}
