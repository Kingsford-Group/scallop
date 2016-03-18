#ifndef __DYNAMIC_GRAPH_H__
#define __DYNAMIC_GRAPH_H__

#include "splice_graph.h"

#include <map>
#include <cassert>

using namespace std;

namespace dynamic_graph
{
	typedef edge_iterator in_edge_iterator;
	typedef edge_iterator out_edge_iterator;

	// basic operations
	int add_vertex(splice_graph &gr);
	PEB add_edge(int s, int t, splice_graph &gr);
	PEB edge(int s, int t, splice_graph &gr);
	int remove_edge(edge_descriptor e, splice_graph &gr);
	int remove_edge(int s, int t, splice_graph &gr);
	int clear_vertex(int x, splice_graph &gr);
	size_t num_vertices(splice_graph &gr);
	size_t num_edges(splice_graph &gr);
	int source(edge_descriptor e, splice_graph &gr);
	int target(edge_descriptor e, splice_graph &gr);
	PEE in_edges(int v, splice_graph &gr);
	PEE out_edges(int v, splice_graph &gr);
	PEE edges(splice_graph &gr);
	int clear(splice_graph &gr);
	int degree(int v, splice_graph &gr);
	int in_degree(int v, splice_graph &gr);
	int out_degree(int v, splice_graph &gr);
	set<int> adjacent_vertices(int v, splice_graph &gr);
	double get_vertex_weight(int v, splice_graph &gr);
	double get_vertex_stddev(int v, splice_graph &gr);
	double get_edge_weight(edge_base *e, splice_graph &gr);
	double get_edge_stddev(edge_base *e, splice_graph &gr);
	int set_vertex_weight(int v, double w, splice_graph &gr);
	int set_vertex_stddev(int v, double w, splice_graph &gr);
	int set_edge_weight(edge_base *e, double w, splice_graph &gr);
	int set_edge_stddev(edge_base *e, double w, splice_graph &gr);

	// read, write, draw and simulate splice graph
	int build_splice_graph(splice_graph &gr, const string &file);
	int write_splice_graph(splice_graph &gr, const string &file);
	int draw_splice_graph(splice_graph &gr, const string &file, double len = 1.2);
	int simulate_splice_graph(splice_graph &gr, int n, int m);

	// get and set properties of splice graph
	int get_edge_weights(splice_graph &gr, MED &med);
	int set_edge_weights(splice_graph &gr, const MED & med);
	int get_vertex_weights(splice_graph &gr, vector<double> &v);
	int set_vertex_weights(splice_graph &gr, const vector<double> &v);
	int get_edge_indices(splice_graph &gr, VE &i2e, MEI &e2i);

	// analysis the structure of splice graph
	int compute_num_paths(splice_graph &gr);
	bool check_nested_splice_graph(splice_graph &gr);
	bool check_fully_connected(splice_graph &gr);
	int bfs_splice_graph(splice_graph &gr, int s, vector<int> &v);
}

#endif
