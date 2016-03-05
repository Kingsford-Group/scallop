#include "decomposer.h"
#include "subsetsum.h"

#include <cstdio>
#include <cmath>
#include <climits>
#include <algorithm>
#include <cassert>

using namespace subsetsum;

decomposer::decomposer(const vector<int> &_s, const vector<int> &_t)
	: s(_s), t(_t)
{
}

int decomposer::solve()
{
	build_subsets();
	//enumerate_paths();
	print();
	return 0;
}

int decomposer::build_subsets()
{
	vector<int> si;
	vector<int> ti;
	for(int i = 0; i < s.size(); i++) si.push_back(i);
	for(int i = 0; i < t.size(); i++) ti.push_back(i);

	PVV p(si, ti);
	vector<PVV> v;
	v.push_back(p);

	int k = 0;
	while(k < v.size())
	{
		vector<int> &xi = v[k].first;
		vector<int> &yi = v[k].second;

		vector<int> x;
		vector<int> y;
		for(int i = 0; i < xi.size(); i++) x.push_back(s[xi[i]]);
		for(int i = 0; i < yi.size(); i++) y.push_back(t[yi[i]]);

		vector<int> subx;
		vector<int> suby;
		double ratio = compute_closest_subsets(x, y, subx, suby);

		if(ratio < 0.99)
		{
			subsets.push_back(v[k]);
		}
		else
		{
			vector<int> s1 = decipher(subx, xi);
			vector<int> t1 = decipher(suby, yi);
			vector<int> s2 = decipher(complement(subx, x.size()), xi);
			vector<int> t2 = decipher(complement(suby, y.size()), yi);
			v.push_back(PVV(s1, t1));
			v.push_back(PVV(s2, t2));
		}

		k++;
	}
	return 0;
}

vector<int> decomposer::complement(const vector<int> &x, int n)
{
	// TODO assume that x is sorted
	vector<int> v;
	int p = 0;
	for(int i = 0; i < x.size(); i++)
	{
		for(int j = p; j != x[i]; j++) v.push_back(j);
		p = x[i] + 1;
	}
	for(int j = p; j < n; j++) v.push_back(j);
	return v;
}

vector<int> decomposer::decipher(const vector<int> &x, const vector<int> &r)
{
	vector<int> v;
	for(int i = 0; i < x.size(); i++)
	{
		v.push_back(r[x[i]]);
	}
	return v;
}

int decomposer::print()
{
	for(int i = 0; i < subsets.size(); i++)
	{
		vector<int> &x = subsets[i].first;
		vector<int> &y = subsets[i].second;
		assert(x.size() >= 1);
		assert(y.size() >= 1);

		printf("subset%d: L = (%d:%d", i, 0, s[x[0]]);
		for(int k = 1; k < x.size(); k++)
		{
			printf(", %d:%d", k, s[x[k]]);
		}
		printf("),  R = (");
		printf("%d:%d", 0, t[y[0]]);
		for(int k = 1; k < y.size(); k++)
		{
			printf(", %d:%d", k, t[y[k]]);
		}
		printf(")\n");
	}

	for(int i = 0; i < vpi.size(); i++)
	{
		printf("decomposition %3d:\n", i);
		for(int j = 0; j < vpi[i].size(); j++)
		{
			PPII &p = vpi[i][j];
			printf(" path %3d: (%3d, %3d) = %4d\n", j, p.first.first, p.first.second, p.second);
		}
	}
	return 0;
}

vector<VPPII> decomposer::enumerate_paths(const vector<int> &x, const vector<int> &y)
{
	vector<VPPII> vp;
	
	if(x.size() == 0 || y.size() == 0) return vp;

	int k = 0;
	int m = INT_MAX;

	for(int i = 1; i < x.size(); i++)
	{
		if(x[i] < m)
		{
			k = i;
			m = x[i];
		}
	}
	for(int i = 0; i < y.size(); i++)
	{
		if(y[i] < m)
		{
			k = i + x.size();
			m = y[i];
		}
	}

	if(k < x.size())
	{
		vector<int> xx;
		for(int i = 0; i < x.size(); i++)
		{
			if(i != k) xx.push_back(x[i]);
		}

		vector<int> yy = y;
		for(int i = 0; i < y.size(); i++)
		{
			assert(yy[i] > x[k]);
			yy[i] -= x[k];

			PPII pi(PI(k, i), x[k]);
			vector<VPPII> v = enumerate_paths(xx, yy);
			for(int j = 0; j < v.size(); j++) v[j].push_back(pi);
			vp.insert(vp.end(), v.begin(), v.end());

			yy[i] += x[k];
		}
	}
	else
	{
		k -= x.size();

		vector<int> yy;
		for(int i = 0; i < y.size(); i++)
		{
			if(i != k) yy.push_back(y[i]);
		}

		vector<int> xx = x;
		for(int i = 0; i < x.size(); i++)
		{
			assert(xx[i] > y[k]);
			xx[i] -= y[k];

			PPII pi(PI(i, k), y[k]);
			vector<VPPII> v = enumerate_paths(xx, yy);
			for(int j = 0; j < v.size(); j++) v[j].push_back(pi);
			vp.insert(vp.end(), v.begin(), v.end());

			xx[i] += y[k];
		}
	}

	return vp;
}

int decomposer::test()
{
	vector<int> s;
	vector<int> t;
	/*
	s.push_back(958);
	s.push_back(189);
	s.push_back(1105);
	t.push_back(443);
	t.push_back(851);
	t.push_back(637);
	t.push_back(319);
	*/

	s.push_back(100);
	s.push_back(20);
	s.push_back(60);
	s.push_back(80);
	t.push_back(100);
	t.push_back(20);
	t.push_back(60);
	t.push_back(80);

	decomposer obj(s, t);
	obj.solve();

	return 0;
}

