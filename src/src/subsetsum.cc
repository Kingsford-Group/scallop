#include "subsetsum.h"
#include "config.h"
#include <cstdio>
#include <cmath>
#include <climits>
#include <algorithm>
#include <cassert>

subsetsum::subsetsum(const vector<int> &v)
	: raw(v)
{
	opt = -1;
	subset.clear();
}

int subsetsum::solve()
{
	if(raw.size() <= 1) return 0;
	init_seeds();
	rescale();
	init_table();
	fill_table();
	optimize();
	backtrace();
	recover();
	//print();
	return 0;
}

int subsetsum::init_seeds()
{
	seeds.clear();
	for(int i = 0; i < raw.size(); i++)
	{
		seeds.push_back(PI(raw[i], i));
	}
	sort(seeds.begin(), seeds.end());
	target = seeds[seeds.size() - 1].first;
	ubound = ceil(target * 1.05);
	seeds.pop_back();
	return 0;
}

int subsetsum::rescale()
{
	// TODO
	return 0;
	int n = seeds.size() * target;
	if(n <= max_dp_table_size) return 0;
	double f = max_dp_table_size * 1.0 / n;
	for(int i = 0; i < seeds.size(); i++)
	{
		seeds[i].first = (int)ceil(seeds[i].first * 1.0 * f);
	}
	target = (int)ceil(target * 1.0 * f);
	ubound = (int)ceil(target * 1.05);
	printf("n = %d, f = %lf, target = %d, ubound = %d\n", n, f, target, ubound);
	return 0;
}

int subsetsum::init_table()
{

	table.resize(seeds.size() + 1);
	btptr.resize(seeds.size() + 1);
	for(int i = 0; i < table.size(); i++)
	{
		table[i].assign(ubound + 1, false);
		btptr[i].assign(ubound + 1, false);
	}

	for(int i = 0; i <= seeds.size(); i++)
	{
		table[i][0] = true;
		btptr[i][0] = false;
	}
	for(int j = 1; j <= ubound; j++)
	{
		table[0][j] = false;
		btptr[0][j] = false;
	}
	return 0;
}

int subsetsum::fill_table()
{
	for(int j = 1; j <= ubound; j++)
	{
		for(int i = 1; i <= seeds.size(); i++)
		{
			int s = seeds[i - 1].first;
			if(j >= s && table[i - 1][j - s] == true)
			{
				table[i][j] = true;
				btptr[i][j] = true;
			}
			else if(table[i - 1][j] == true)
			{
				table[i][j] = true;
				btptr[i][j] = false;
			}
			else
			{
				table[i][j] = false;
				btptr[i][j] = false;
			}
		}
	}
	return 0;
}

int subsetsum::optimize()
{
	int s = seeds.size();
	int ku = INT_MAX;
	for(int i = target + 1; i <= ubound; i++)
	{
		if(table[s][i] == true) 
		{
			ku = i;
			break;
		}
	}
	int kl = -1;
	for(int i = target; i >= 1; i--)
	{
		if(table[s][i] == true)
		{
			kl = i;
			break;
		}
	}
	assert(kl != -1);

	opt = kl;
	if(ku - target < target - kl) opt = ku;
	return 0;
}

int subsetsum::backtrace()
{
	int s = seeds.size();
	int x = opt;
	while(x >= 1 && s >= 1)
	{
		assert(table[s][x] == true);
		bool b = btptr[s][x];
		if(b == true)
		{
			x = x - seeds[s - 1].first;
			subset.push_back(s - 1);
		}
		s = s - 1;
	}
	return 0;
}

int subsetsum::recover()
{
	vector<int> v;
	int sum = 0;
	for(int i = 0; i < subset.size(); i++)
	{
		int k = subset[i];
		int r = seeds[k].second;
		v.push_back(r);
		sum += raw[r];
	}
	opt = sum;
	subset = v;
	return 0;
}

int subsetsum::print()
{
	printf(" subsetsum = (%d", raw[0]);
	for(int i = 1; i < raw.size(); i++) printf(", %d", raw[i]);
	printf("), opt = %d, subset = (%d:%d", opt, subset[0], raw[subset[0]]);
	for(int i = 1; i < subset.size(); i++) printf(", %d:%d", subset[i], raw[subset[i]]);
	printf(")\n");
	return 0;
}

int subsetsum::test()
{
	vector<int> v;
	v.push_back(100);
	v.push_back(200);
	v.push_back(1);
	v.push_back(3);
	v.push_back(150);

	subsetsum sss(v);
	sss.solve();
	printf("opt = %d, subset = (%d", sss.opt, sss.subset[0]);
	for(int i = 1; i < sss.subset.size(); i++) printf(", %d", sss.subset[i]);
	printf(")\n");
	return 0;
}


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

