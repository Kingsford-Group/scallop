#ifndef __SPLICE_GRAPH_H__
#define __SPLICE_GRAPH_H__

#include "directed_graph.h"

#include <map>
#include <cassert>

#define SMIN 0.0000000001

using namespace std;

class splice_graph : public directed_graph
{
public:
	splice_graph();
	splice_graph(const splice_graph &gr);
	virtual ~splice_graph();

private:
	vector<string> vstr;
	vector<double> vwrt;
	vector<double> vdev;
	MED ewrt;
	MED edev;

public:
	// get and set properties
	string get_vertex_string(int v) const;
	double get_vertex_weight(int v) const;
	double get_vertex_stddev(int v) const;
	double get_edge_weight(edge_base *e) const;
	double get_edge_stddev(edge_base *e) const;

	int set_vertex_string(int v, string s);
	int set_vertex_weight(int v, double w);
	int set_vertex_stddev(int v, double w);
	int set_edge_weight(edge_base *e, double w);
	int set_edge_stddev(edge_base *e, double w);

	MED get_edge_weights() const;
	vector<double> get_vertex_weights() const;
	int set_edge_weights(const MED & med);
	int set_vertex_weights(const vector<double> &v);

	// modify the splice_graph
	int clear();

	// read, write, and simulate splice graph
	int build(const string &file);
	int write(const string &file) const;
	int simulate(int n, int m);

	// analysis the structure of splice graph
	int compute_num_paths();
	int compute_decomp_paths();
	bool check_fully_connected();

	// algorithms with weight contraints
	int bfs_w(int s, double w, vector<int> &v, VE &b);
	int compute_shortest_path_w(int s, int t, double w);
	int compute_shortest_path_w(int s, int t, double w, VE &p);
	double compute_maximum_path_w(VE &p);
	double compute_minimum_weight(const VE &p);

	// determine optimal path
	bool compute_optimal_path(VE &p);

	// rounding all weights to integers
	int round_weights();

	// remove edges with 0 weight
	int remove_empty_edges();

	// draw
	int draw(const string &file);
	int draw(const string &file, const MIS &mis, const MES &mes, double len);
};

#endif
