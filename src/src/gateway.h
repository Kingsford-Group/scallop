#ifndef __GATEWAY_H__
#define __GATEWAY_H__

#include <vector>

using namespace std;

typedef pair<int, int> PI;

class gateway
{
public:
	vector<PI> routes;
	vector<double> counts;

public:
	int add_route(const PI &p, double c);
	int remove_route(const PI &p);
	int replace_in_edge(int ex, int ey);
	int replace_out_edge(int ex, int ey);
	int split_in_edge(int ex, int ey, double r);
	int split_out_edge(int ex, int ey, double r);
	int remove_in_edges(const vector<int> &v);
	int remove_out_edges(const vector<int> &v);
	int remove_in_edge(int x);
	int remove_out_edge(int x);
	int print(int index) const;
	double total_counts() const;
};

#endif
