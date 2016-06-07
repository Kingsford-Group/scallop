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
	optimize();

	for(int i = 0; i < opts.size(); i++)
	{
		vector<int> subset;
		backtrace(opts[i], subset);
		recover(opts[i], subset);
		subsets.push_back(subset);
	}
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
			if(j >= s && table[i - 1][j - s] == 2)
			{
				table[i][j] = -1;
			}
			else if(table[i - 1][j] == 2)
			{
				table[i][j] = -1;
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

int subsetsum2::optimize()
{
	int s = seeds.size();
	bool b = true;
	int d = 0;
	while(true)
	{
		int k = target + (b ? 0 + d : 0 - d);
		if(k <= 0 || k > ubound) break;
		//printf("d = %d, target = %d, ubound = %d, table[%d][%d] = %d\n", d, target, ubound, s, k, table[s][k] ? 1 : 0);
		if(table[s][k] == true)
		{
			opts.push_back(k);
			if(opts.size() >= max_num_subsetsum_solutions) break;
		}
		if(b == true) d++;
		b = (b ? false : true);
	}
	//assert(opts.size() >= 1);
	return 0;
}

int subsetsum2::backtrace(int opt, vector<int> &subset)
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

int subsetsum2::recover(int &opt, vector<int> &subset)
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

