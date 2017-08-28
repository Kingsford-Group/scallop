/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

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
	directed_graph(const directed_graph &gr);
	directed_graph& operator=(const directed_graph &gr);
	virtual ~directed_graph();

public:
	// modify the graph
	virtual edge_descriptor add_edge(int s, int t);
	virtual int remove_edge(edge_descriptor e);
	virtual int remove_edge(int s, int t);
	virtual int move_edge(edge_base *e, int x, int y);
	virtual int exchange(int x, int y, int z);
	virtual int rotate(int x, int y);

	// access functions
	virtual int in_degree(int v) const;
	virtual int out_degree(int v) const;
	virtual PEEI in_edges(int v);
	virtual PEEI out_edges(int v);

	// algorithms
	virtual int bfs_reverse(int t, vector<int> &v);
	virtual int bfs_reverse(int t, vector<int> &v, vector<int> &b);
	virtual int bfs_reverse(int t, set<edge_descriptor> &ss);
	virtual bool bfs_reverse(const vector<int> &vt, int s, const set<edge_descriptor> &fb);
	virtual bool compute_shortest_path(int x, int y, vector<int> &p);
	virtual bool compute_shortest_path(edge_descriptor ex, edge_descriptor ey, vector<int> &p);
	virtual bool check_path(int x, int y);
	virtual bool check_path(edge_descriptor ex, edge_descriptor ey);
	virtual bool intersect(edge_descriptor ex, edge_descriptor ey);
	virtual vector<int> topological_sort();
	virtual vector<int> topological_sort_reverse();
	virtual vector<int> topological_sort0();
	virtual int compute_in_partner(int x);
	virtual int compute_out_partner(int x);
	virtual int compute_in_equivalent_vertex(int x);
	virtual int compute_out_equivalent_vertex(int x);
	virtual int check_nest(int x, int r, set<edge_descriptor> &vv);
	virtual int check_nest(int x, int r, set<edge_descriptor> &vv, const vector<int> &tpo);
	virtual int check_nest(int x, int r, const vector<int> &tpo);

	// draw
	int draw(const string &file, const MIS &mis, const MES &mes, double len);
	int draw(const string &file, const MIS &mis, const MES &mes, double len, const vector<int> &topo);
};

#endif
