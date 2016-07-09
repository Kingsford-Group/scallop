#include "subsetsum1.h"
#include "config.h"
#include <cstdio>
#include <cmath>
#include <climits>
#include <algorithm>
#include <cassert>

subsetsum1::subsetsum1(const vector<int> &s, const vector<int> &t)
{
	for(int i = 0; i < s.size(); i++)
	{
		source.push_back(PI(s[i], i));
	}
	for(int i = 0; i < t.size(); i++)
	{
		target.push_back(PI(t[i], i));
	}
	eqns.clear();
}

int subsetsum1::solve()
{
	init();
	fill();
	optimize();
	return 0;
}

int subsetsum1::init()
{
	sort(source.begin(), source.end());
	sort(target.begin(), target.end());
	ubound = target[target.size() - 1].first;

	table1.resize(source.size() + 1);
	table2.resize(source.size() + 1);
	for(int i = 0; i < table1.size(); i++)
	{
		table1[i].assign(ubound + 1, false);
		table2[i].assign(ubound + 1, -1);
	}

	for(int i = 0; i <= source.size(); i++)
	{
		table1[i][0] = false;
		table2[i][0] = 0;
	}

	for(int j = 1; j <= ubound; j++)
	{
		table1[0][j] = false;
		table2[0][j] = -1;
	}
	return 0;
}

int subsetsum1::fill()
{
	for(int j = 1; j <= ubound; j++)
	{
		for(int i = 1; i <= source.size(); i++)
		{
			int s = source[i - 1].first;
			if(j >= s && table2[i - 1][j - s] >= 0)
			{
				table1[i][j] = true;
				table2[i][j] = i;
			}
			
			if(table2[i][j] == -1 && table2[i - 1][j] >= 0)
			{
				table2[i][j] = table2[i - 1][j];
			}
		}
	}
	return 0;
}

int subsetsum1::optimize()
{
	eqns.clear();
	for(int i = 0; i < target.size(); i++)
	{
		equation eqn(0);
		backtrace(i, eqn);
		if(eqn.s.size() == 0) continue;
		eqns.push_back(eqn);
	}
	return 0;
}

int subsetsum1::backtrace(int ti, equation &eqn)
{
	int t = target[ti].first;
	int n = source.size();
	if(table2[n][t] == -1) return 0;

	eqn.s.push_back(target[ti].second);

	int x = t;
	int s = table2[n][t];
	while(x >= 1 && s >= 1)
	{
		assert(table2[s][x] >= 0);
		eqn.t.push_back(source[s - 1].second);

		x -= source[s - 1].first;
		s = table2[s - 1][x];
	}
	return 0;
}

int subsetsum1::print()
{
	printf("table 1\n");
	printf("   ");
	for(int i = 0; i < table1[0].size(); i++) printf("%3d", i);
	printf("\n");

	for(int i = 0; i < table1.size(); i++)
	{
		printf("%3d", i);
		for(int j = 0; j < table1[i].size(); j++)
		{
			printf("%3d", table1[i][j] ? 1 : 0);
		}
		printf("\n");
	}

	printf("table 2\n");
	printf("   ");
	for(int i = 0; i < table2[0].size(); i++) printf("%3d", i);
	printf("\n");

	for(int i = 0; i < table2.size(); i++)
	{
		printf("%3d", i);
		for(int j = 0; j < table2[i].size(); j++)
		{
			printf("%3d", table2[i][j]);
		}
		printf("\n");
	}

	for(int i = 0; i < eqns.size(); i++)
	{
		eqns[i].print(i);
	}
	return 0;
}

int subsetsum1::test()
{
	vector<int> v;
	v.push_back(5);
	v.push_back(6);
	v.push_back(8);
	v.push_back(3);
	v.push_back(10);

	vector<int> t;
	t.push_back(20);
	t.push_back(8);
	t.push_back(11);
	t.push_back(18);

	subsetsum1 sss(v, t);
	sss.solve();
	sss.print();


	return 0;
}
