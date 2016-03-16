#include "subsetsum.h"
#include <cstdio>
#include <cmath>
#include <climits>
#include <algorithm>
#include <cassert>

subsetsum::subsetsum(const vector<int> &v)
	: seeds(v)
{}

double subsetsum::compute_closest_subsets(const vector<int> &s, const vector<int> &t, vector<int> &subs, vector<int> &subt)
{
	if(s.size() <= 1 || t.size() <= 1) return 0;

	vector<int> ss;
	vector<int> tt;
	vector<int> sf;
	vector<int> tf;
	vector<int> sb;
	vector<int> tb;
	enumerate_subsets(s, ss, sf, sb);
	enumerate_subsets(t, tt, tf, tb);

	int ssi;
	int tti;
	int dist = compute_closest_pair(ssi, tti, ss, tt);

	recover_subset(subs, ssi, sf, sb);
	recover_subset(subt, tti, tf, tb);

	sort(subs.begin(), subs.end());
	sort(subt.begin(), subt.end());

	compute_subset_ratio(s, t, subs, subt);

	return 0;
}

int subsetsum::enumerate_subsets(const vector<int> &x, vector<int> &xx)
{
	vector<int> xf;
	vector<int> xb;
	enumerate_subsets(x, xx, xf, xb);
	return 0;
}

int subsetsum::enumerate_subsets(const vector<int> &x, vector<int> &xx, vector<int> &xf, vector<int> &xb)
{
	xx = x;
	int start = 0;
	xf.clear();
	xb.clear();
	for(int i = 0; i < x.size(); i++) xf.push_back(i + 1);
	for(int i = 0; i < x.size(); i++) xb.push_back(-1);
	for(int k = 0; k < x.size() - 2; k++) augment(x, xx, xf, xb, start);
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

int subsetsum::compute_closest_pair(int &ssi, int &tti, const vector<int> &ss, const vector<int> &tt)
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
	int dd = INT_MAX;

	while(si < sss.size() && ti < ttt.size())
	{
		int d = (int)(fabs(sss[si].first - ttt[ti].first));
		if(d < dd)
		{
			dd = d;
			ssi = sss[si].second;
			tti = ttt[ti].second;
		}

		if(dd == 0) break;
		if(sss[si].first < ttt[ti].first) si++;
		else ti++;
	}
	return dd;
}

int subsetsum::recover_subset(vector<int> &sub, int xxi, const vector<int> &xf, const vector<int> &xb)
{
	sub.clear();
	int x = xxi;
	while(true)
	{
		assert(x >= 0 && x < xf.size());
		assert(xf[x] >= 1);
		sub.push_back(xf[x] - 1);
		if(xb[x] == -1) break;
		x = xb[x];
	}
	return 0;
}

double subsetsum::compute_subset_ratio(const vector<int> &s, const vector<int> &t, const vector<int> &subs, const vector<int> &subt)
{
	int ss = 0, tt = 0;
	for(int i = 0; i < s.size(); i++) ss += s[i];
	for(int i = 0; i < t.size(); i++) tt += t[i];
	double x0 = ss < tt ? ss : tt;

	int ss1 = 0, tt1 = 0;
	for(int i = 0; i < subs.size(); i++) ss1 += s[subs[i]];
	for(int i = 0; i < subt.size(); i++) tt1 += t[subt[i]];

	int ss2 = ss - ss1;
	int tt2 = tt - tt1;

	assert(ss1 > 0 && ss2 > 0);
	assert(tt1 > 0 && tt2 > 0);

	double x1 = ss1 < tt1 ? ss1 : tt1;
	double x2 = ss2 < tt2 ? ss2 : tt2;

	return (x1 + x2) / x0;
}

