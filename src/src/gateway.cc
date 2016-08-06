#include "gateway.h"
#include "util.h"
#include <cstdio>
#include <algorithm>
#include <set>

int gateway::print(int index) const
{
	printf("gateway %d, #routes = %lu, total-counts = %.1lf\n", index, routes.size(), total_counts());
	for(int i = 0; i < routes.size(); i++)
	{
		printf(" route %d (%d, %d), count = %.1lf\n", i, routes[i].first, routes[i].second, counts[i]);
	}
	return 0;
}

double gateway::total_counts() const
{
	double s = 0;
	for(int i = 0; i < counts.size(); i++)
	{
		s += counts[i];
	}
	return s;
}

int gateway::add_route(const PI &p, double c)
{
	for(int i = 0; i < routes.size(); i++)
	{
		if(routes[i] == p)
		{
			counts[i] += c;
			return 0;
		}
	}

	routes.push_back(p);
	counts.push_back(c);
	return 0;
}

int gateway::remove_route(const PI &p)
{
	for(int i = 0; i < routes.size(); i++)
	{
		if(routes[i] == p)
		{
			routes.erase(routes.begin() + i);
			counts.erase(counts.begin() + i);
			return 0;
		}
	}
	return 0;
}

int gateway::replace_in_edge(int ex, int ey)
{
	for(int i = 0; i < routes.size(); i++)
	{
		if(routes[i].first == ex) routes[i].first = ey;
	}
	return 0;
}

int gateway::replace_out_edge(int ex, int ey)
{
	for(int i = 0; i < routes.size(); i++)
	{
		if(routes[i].second == ex) routes[i].second = ey;
	}
	return 0;
}

int gateway::split_in_edge(int ex, int ey, double r)
{
	assert(r >= 0 && r <= 1.0);
	if(ex == ey) return 0;
	int n = routes.size();
	for(int i = 0; i < n; i++)
	{
		if(routes[i].first == ex)
		{
			add_route(PI(ey, routes[i].second), (1.0 - r) * counts[i]);
			counts[i] *= r;
		}
	}
	return 0;
}

int gateway::split_out_edge(int ex, int ey, double r)
{
	assert(r >= 0 && r <= 1.0);
	if(ex == ey) return 0;
	int n = routes.size();
	for(int i = 0; i < n; i++)
	{
		if(routes[i].second == ex)
		{
			add_route(PI(routes[i].first, ey), (1.0 - r) * counts[i]);
			counts[i] *= r;
		}
	}
	return 0;
}

int gateway::remove_in_edges(const vector<int> &v)
{
	set<int> s(v.begin(), v.end());

	vector<PI> vv;
	vector<double> cc;
	for(int i = 0; i < routes.size(); i++)
	{
		int x = routes[i].first;
		if(s.find(x) != s.end()) continue;
		vv.push_back(routes[i]);
		cc.push_back(counts[i]);
	}

	routes = vv;
	counts = cc;

	return 0;
}

int gateway::remove_out_edges(const vector<int> &v)
{
	set<int> s(v.begin(), v.end());

	vector<PI> vv;
	vector<double> cc;
	for(int i = 0; i < routes.size(); i++)
	{
		int x = routes[i].second;
		if(s.find(x) != s.end()) continue;
		vv.push_back(routes[i]);
		cc.push_back(counts[i]);
	}

	routes = vv;
	counts = cc;

	return 0;
}

int gateway::remove_in_edge(int x)
{
	vector<int> v;
	v.push_back(x);
	return remove_in_edges(v);
}

int gateway::remove_out_edge(int x)
{
	vector<int> v;
	v.push_back(x);
	return remove_out_edges(v);
}
