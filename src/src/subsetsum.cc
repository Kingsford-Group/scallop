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
	opts.clear();
	subsets.clear();
}

int subsetsum::solve()
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

int subsetsum::init_seeds()
{
	seeds.clear();
	for(int i = 0; i < raw.size(); i++)
	{
		seeds.push_back(PI(raw[i], i));
	}
	sort(seeds.begin(), seeds.end());
	target = seeds[seeds.size() - 1].first;
	ubound = ceil(target * 2.0);
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
	bool b = true;
	int d = 0;
	while(true)
	{
		int k = target + (b ? 0 + d : 0 - d);
		if(k <= 1 || k >= ubound) break;
		if(table[s][k] == true)
		{
			opts.push_back(k);
			if(opts.size() >= max_num_subsetsum_solutions) break;
		}
		b = (b ? false : true);
		d++;
	}
	assert(opts.size() >= 1);
	return 0;
}

int subsetsum::backtrace(int opt, vector<int> &subset)
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

int subsetsum::recover(int &opt, vector<int> &subset)
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

