#include "subsetsum2.h"
#include "config.h"
#include <cstdio>
#include <cmath>
#include <climits>
#include <algorithm>
#include <cassert>

subsetsum2::subsetsum2(const vector<int> &v)
	: raw(v)
{
	eqns.clear();
}

int subsetsum2::solve()
{
	if(raw.size() <= 1) return 0;
	init_seeds();
	rescale();
	init_table();
	fill_table();
	optimize1();
	optimize2();
	recover();
	return 0;
}

int subsetsum2::init_seeds()
{
	seeds.clear();
	int sum = 0;
	for(int i = 0; i < raw.size(); i++)
	{
		sum += raw[i];
		seeds.push_back(PI(raw[i], i));
	}
	sort(seeds.begin(), seeds.end());
	ubound = (sum + 1) / 2;
	return 0;
}

int subsetsum2::rescale()
{
	// TODO
	return 0;
	int n = seeds.size() * ubound;
	if(n <= max_dp_table_size) return 0;
	double f = max_dp_table_size * 1.0 / n;
	int sum = 0;
	for(int i = 0; i < seeds.size(); i++)
	{
		seeds[i].first = (int)ceil(seeds[i].first * 1.0 * f);
		sum += seeds[i].first;
	}
	ubound = (sum + 1) / 2 * 1.05;
	return 0;
}

int subsetsum2::init_table()
{
	table.resize(seeds.size() + 1);
	for(int i = 0; i < table.size(); i++)
	{
		table[i].assign(ubound + 1, -1);
	}

	for(int i = 0; i <= seeds.size(); i++)
	{
		table[i][0] = 0;
	}
	for(int j = 1; j <= ubound; j++)
	{
		table[0][j] = -1;
	}
	return 0;
}

int subsetsum2::fill_table()
{
	for(int j = 1; j <= ubound; j++)
	{
		for(int i = 1; i <= seeds.size(); i++)
		{
			int s = seeds[i - 1].first;
			if(j >= s && table[i - 1][j - s] >= 2)
			{
				table[i][j] = 4;
			}
			else if(table[i - 1][j] >= 2)
			{
				table[i][j] = 3;
			}
			else if(j >= s && table[i - 1][j - s] >= 0 && table[i - 1][j] >= 0)
			{
				table[i][j] = 2;
			}
			else if(j >= s && table[i - 1][j - s] >= 0)
			{
				table[i][j] = 1;
			}
			else if(table[i - 1][j] >= 0)
			{
				table[i][j] = 0;
			}
			else 
			{
				table[i][j] = -1;
			}
		}
	}
	return 0;
}

int subsetsum2::optimize1()
{
	// perfect equation (error = 0)
	for(int j = 1; j < ubound; j++)
	{
		int i = 0;
		for(int k = 1; k <= seeds.size(); k++)
		{
			if(table[k][j] == 2)
			{
				i = k;
				break;
			}
		}

		if(i == 0) continue;
		//printf("verify (%d, %d)\n", i, j);
		if(verify(i, j) == false) continue;

		equation eqn(0);
		backtrace(i - 1, j, eqn.s);
		vector<int> v;
		backtrace(i - 1, j - seeds[i - 1].first, v);
		eqn.t.push_back(i - 1);
		eqn.t.insert(eqn.t.end(), v.begin(), v.end());
		eqns.push_back(eqn);
	}
	return 0;
}

bool subsetsum2::verify(int s, int x)
{
	for(int i = 0; i < s; i++)
	{
		int xx = x - seeds[i].first;
		if(table[i + 1][xx] >= 2) return false;
	}
	return true;
}

int subsetsum2::optimize2()
{
	// TODO
	// non-perfect equation (error > 0)
	vector<int> v;
	int s = seeds.size();
	for(int x = 1; x <= ubound; x++)
	{
		if(table[s][x] < 0) continue;
		if(table[s][x] >= 2) continue;
		v.push_back(x);
	}

	if(v.size() == 0) return 0;

	vector<int> s1;
	vector<int> s2;
	backtrace(s, v[0], s1);
	for(int i = 1; i < v.size(); i++)
	{
		backtrace(s, v[i], s2);
		double r = (v[i] - v[i - 1]) * 1.0 / v[i - 1];
		if(r >= max_equation_error_ratio)
		{
			s1 = s2;
			continue;
		}
		equation eqn(s1, s2, v[i] - v[i - 1]);
		eqns.push_back(eqn);
	}
	return 0;
}

int subsetsum2::backtrace(int s, int x, vector<int> &subset)
{
	subset.clear();
	while(x >= 1 && s >= 1)
	{
		if(table[s][x] == 0)
		{
			s = s - 1;
		}
		else if(table[s][x] == 1)
		{
			s = s - 1;
			x = x - seeds[s].first;
			subset.push_back(s);
		}
		else
		{
			assert(false);
		}
	}
	return 0;
}

int subsetsum2::recover()
{
	for(int i = 0; i < eqns.size(); i++)
	{
		int s1 = recover(eqns[i].s);
		int s2 = recover(eqns[i].t);
		eqns[i].e = (int)fabs(s1 - s2);
	}
	return 0;
}

int subsetsum2::recover(vector<int> &subset)
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
	return sum;
}

int subsetsum2::test()
{
	vector<int> v;
	v.push_back(1); // 0
	v.push_back(2); // 1
	v.push_back(3); // 2
	v.push_back(4); // 3
	v.push_back(5); // 4
	v.push_back(6); // 5
	v.push_back(7); // 6
	v.push_back(8); // 7
	v.push_back(9); // 8

	subsetsum2 sss(v);
	sss.solve();

	for(int i = 0; i < sss.eqns.size(); i++)
	{
		sss.eqns[i].print(i);
	}

	return 0;
}
int subsetsum2::print()
{
	int s = seeds.size();
	for(int i = 1; i <= ubound; i++)
	{
		printf("table[%d] = %d\n", i, table[s][i]);
	}
	printf("\n");
	return 0;
}
