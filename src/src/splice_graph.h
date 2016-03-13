#ifndef __SPLICE_GRAPH_H__
#define __SPLICE_GRAPH_H__

#include "graph_base.h"

#include <map>
#include <cassert>

using namespace std;

class splice_graph : public graph_base
{
public:
	splice_graph();
	splice_graph(const splice_graph &gr);
	virtual ~splice_graph();

private:
	vector<double> vwrt;
	vector<double> vdev;
	MED ewrt;
	MED edev;

public:
	// get and set properties
	double get_vertex_weight(int v) const;
	double get_vertex_stddev(int v) const;
	double get_edge_weight(edge_base *e) const;
	double get_edge_stddev(edge_base *e) const;
	int set_vertex_weight(int v, double w);
	int set_vertex_stddev(int v, double w);
	int set_edge_weight(edge_base *e, double w);
	int set_edge_stddev(edge_base *e, double w);

	MED get_edge_weights() const;
	vector<double> get_vertex_weights() const;
	int set_edge_weights(const MED & med);
	int set_vertex_weights(const vector<double> &v);
	int get_edge_indices(VE &i2e, MEI &e2i) const;

	// modify the splice_graph
	int clear();

	// read, write, draw and simulate splice graph
	int build(const string &file);
	int write(const string &file) const;
	int draw(const string &file, double len = 1.2) const;
	int simulate(int n, int m);

	// analysis the structure of splice graph
	int compute_num_paths() const;
	bool check_fully_connected() const;
	bool check_nested() const;
};

#endif
