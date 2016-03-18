#ifndef __DIRECTED_GRAPH_H__
#define __DIRECTED_GRAPH_H__

#include <vector>
#include <map>

#include "graph_base.h"

using namespace std;

class directed_graph : public graph_base
{
public:
	directed_graph();
	directed_graph(const graph_base &gr);
	virtual ~directed_graph();

public:
	// modify the graph
	virtual edge_descriptor add_edge(int s, int t);
	virtual int remove_edge(edge_descriptor e);
	virtual int remove_edge(int s, int t);
	virtual int move_edge(edge_base *e, int x, int y);
	virtual int exchange(int x, int y, int z);

	// access functions
	virtual int in_degree(int v) const;
	virtual int out_degree(int v) const;
	virtual PEE in_edges(int v);
	virtual PEE out_edges(int v);

	// algorithms
	virtual int bfs_reverse(int t, vector<int> &v);
	virtual int bfs_reverse(int t, vector<int> &v, vector<int> &b);
	virtual bool compute_shortest_path(int x, int y, vector<int> &p);
	virtual bool compute_shortest_path(edge_descriptor ex, edge_descriptor ey, vector<int> &p);
	virtual bool check_path(int x, int y);
	virtual bool check_path(edge_descriptor ex, edge_descriptor ey);
	virtual bool intersect(edge_descriptor ex, edge_descriptor ey);
	virtual vector<int> topological_sort();
	virtual int compute_in_partner(int x);
	virtual int compute_out_partner(int x);

	// print and draw
	int draw(const string &file, const MIS &mis, const MES &mes, double len);
};

#endif
