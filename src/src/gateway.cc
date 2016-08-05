#include "gateway.h"
#include "util.h"
#include <cstdio>
#include <algorithm>

int gateway::print(int index) const
{
	printf("gateway %d, #routes = %lu, total-counts = %d\n", index, routes.size(), total_counts());
	for(int i = 0; i < routes.size(); i++)
	{
		printf(" route %d (%d, %d), count = %d\n", i, routes[i].first, routes[i].second, counts[i]);
	}
	return 0;
}

int gateway::total_counts() const
{
	int s = 0;
	for(int i = 0; i < counts.size(); i++)
	{
		s += counts[i];
	}
	return s;
}

int gateway::add_route(const PI &p, int c)
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
