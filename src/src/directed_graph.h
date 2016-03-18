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
	virtual int remove_edge(edge_descriptor e);
	virtual int remove_edge(int s, int t);
	virtual int move_edge(edge_base *e, int x, int y);
	virtual int exchange(int x, int y, int z);

	// access functions
	virtual PEB edge(int s, int t) const;
	virtual int in_degree(int v) const;
	virtual int out_degree(int v) const;
	virtual PEE edges() const;
	virtual vector<edge_descriptor> edges(int x, int y) const;
	virtual PEE in_edges(int v) const;
	virtual PEE out_edges(int v) const;
	virtual set<int> adjacent_vertices(int v) const;

	// algorithms
	virtual int bfs_reverse(int t, vector<int> &v) const;
	virtual int bfs_reverse(int t, vector<int> &v, vector<int> &b) const;
	virtual bool compute_shortest_path(int x, int y, vector<int> &p) const;
	virtual bool compute_shortest_path(edge_descriptor ex, edge_descriptor ey, vector<int> &p) const;
	virtual bool check_path(int x, int y) const;
	virtual bool check_path(edge_descriptor ex, edge_descriptor ey) const;
	virtual bool intersect(edge_descriptor ex, edge_descriptor ey) const;
	virtual vector<int> topological_sort() const;
	virtual int compute_in_partner(int x) const;
	virtual int compute_out_partner(int x) const;

	// print and draw
	int draw(const string &file, const MIS &mis, const MES &mes, double len) const;
	int print() const;
};

#endif
