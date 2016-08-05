#ifndef __GATEWAY_H__
#define __GATEWAY_H__

#include <vector>

using namespace std;

typedef pair<int, int> PI;

class gateway
{
public:
	vector<PI> routes;
	vector<int> counts;

public:
	int add_route(const PI &p, int c);
	int replace_in_edge(int ex, int ey);
	int replace_out_edge(int ex, int ey);
	int print(int index) const;
	int total_counts() const;
};

#endif
