#ifndef __DYNAMIC_GRAPH_H__
#define __DYNAMIC_GRAPH_H__

#include "graph_base.h"

#include <map>
#include <cassert>

using namespace std;

namespace dynamic_graph
{
	typedef edge_descriptor_b edge_descriptor;
	typedef edge_iterator_b edge_iterator;
	typedef edge_iterator_b in_edge_iterator;
	typedef edge_iterator_b out_edge_iterator;
	typedef adj_iterator_b adj_iterator;
	typedef PAA_b PAA;
	typedef PEE_b PEE;
	typedef PEB_b PEB;

	typedef map<edge_descriptor, bool> MEB;
	typedef map<edge_descriptor, double> MED;
	typedef pair<edge_descriptor, double> PED;
	typedef map<edge_descriptor, int> MEI;
	typedef pair<edge_descriptor, int> PEI;
	typedef vector<edge_descriptor> VE;


	class splice_graph : public graph_b
	{
	public:
		splice_graph();
		virtual ~splice_graph();

	public:
		vector<double> vwrt;
		vector<double> vdev;
		map<edge_b*, double> ewrt;
		map<edge_b*, double> edev;

	public:
		double get_vertex_weight(int v) const;
		double get_vertex_stddev(int v) const;
		double get_edge_weight(edge_b *e) const;
		double get_edge_stddev(edge_b *e) const;
		int set_vertex_weight(int v, double w);
		int set_vertex_stddev(int v, double w);
		int set_edge_weight(edge_b *e, double w);
		int set_edge_stddev(edge_b *e, double w);
	};

	// wrappers
	int add_vertex(splice_graph &gr);
	PEB add_edge(int s, int t, splice_graph &gr);
	PEB edge(int s, int t, const splice_graph &gr);
	int remove_edge(edge_descriptor e, splice_graph &gr);
	int remove_edge(int s, int t, splice_graph &gr);
	int clear_vertex(int x, splice_graph &gr);
	size_t num_vertices(const splice_graph &gr);
	size_t num_edges(const splice_graph &gr);
	int source(edge_descriptor e, const splice_graph &gr);
	int target(edge_descriptor e, const splice_graph &gr);
	PEE in_edges(int v, const splice_graph &gr);
	PEE out_edges(int v, const splice_graph &gr);
	PEE edges(const splice_graph &gr);
	int clear(splice_graph &gr);
	int degree(int v, const splice_graph &gr);
	int in_degree(int v, const splice_graph &gr);
	int out_degree(int v, const splice_graph &gr);
	PAA adjacent_vertices(int v, const splice_graph &gr);
	double get_vertex_weight(int v, const splice_graph &gr);
	double get_vertex_stddev(int v, const splice_graph &gr);
	double get_edge_weight(edge_b *e, const splice_graph &gr);
	double get_edge_stddev(edge_b *e, const splice_graph &gr);
	int set_vertex_weight(int v, double w, splice_graph &gr);
	int set_vertex_stddev(int v, double w, splice_graph &gr);
	int set_edge_weight(edge_b *e, double w, splice_graph &gr);
	int set_edge_stddev(edge_b *e, double w, splice_graph &gr);

	// read, write, draw and simulate splice graph
	int build_splice_graph(splice_graph &gr, const string &file);
	int write_splice_graph(const splice_graph &gr, const string &file);
	int draw_splice_graph(const splice_graph &gr, const string &file, double len = 1.2);
	int simulate_splice_graph(splice_graph &gr, int n, int m);

	// get and set properties of splice graph
	int get_edge_weights(const splice_graph &gr, MED &med);
	int set_edge_weights(splice_graph &gr, const MED & med);
	int get_vertex_weights(const splice_graph &gr, vector<double> &v);
	int set_vertex_weights(splice_graph &gr, const vector<double> &v);
	int get_edge_indices(const splice_graph &gr, VE &i2e, MEI &e2i);

	// analysis the structure of splice graph
	int compute_num_paths(const splice_graph &gr);
	bool check_nested_splice_graph(const splice_graph &gr);
	bool check_directed_path(const splice_graph &gr, int s, int t);
	bool check_fully_connected(const splice_graph &gr);
	bool check_fully_reachable_from_start(const splice_graph &gr);
	bool check_fully_reachable_to_end(const splice_graph &gr);
	int bfs_distance(const splice_graph &gr, int s, vector<int> &v);

	// tests
	int test_bfs_distance();
	int test_remove_edge();
}

#endif
