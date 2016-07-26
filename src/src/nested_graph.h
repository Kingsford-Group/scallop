#ifndef __NESTED_GRAPH_H__
#define __NESTED_GRAPH_H__

#include "directed_graph.h"

#include <map>
#include <cassert>
#include <string>

using namespace std;

typedef pair<int, int> PI;
typedef pair<int, edge_descriptor> PIE;

class nested_graph : public directed_graph
{
public:
	nested_graph();
	nested_graph(directed_graph &gr);
	virtual ~nested_graph();

public:
	vector<int> order;		// partial order of vertex
	vector<PEE> partners;	// in/out partners for each vertex
	MEE parents;			// parent for each edge

public:
	int clear();
	int build(directed_graph &gr);
	int draw(const string &file);
	bool link(int xs, int xt, int ys, int yt);
	bool link(int xs, int xt, int ys, int yt, vector<PI> &xp, vector<PI> &yp);

private:
	int init(directed_graph &gr);
	int build_nests0(directed_graph &gr);
	int build_nests(directed_graph &gr);
	bool verify_minimal_nest(directed_graph &gr, const vector<int> &tpo, vector< set<int> > &vs, vector< set<int> > &vt, int i, int j);
	int build_partial_order();
	int build_partners(directed_graph &gr);
	int build_parents();
	int build_parents(int x);

	edge_descriptor get_left_sibling(edge_descriptor e);
	edge_descriptor get_right_sibling(edge_descriptor e);

	int move_left(edge_descriptor e);
	int move_right(edge_descriptor e);

	bool dock_left(edge_descriptor e, vector<int> &r);
	bool dock_right(edge_descriptor e, vector<int> &r);

	edge_descriptor compute_lca(edge_descriptor x, edge_descriptor y);
	int dock(edge_descriptor e, edge_descriptor p, vector<PI> &v);

	PEB add_extra_edge(int s, int t);
	int remove_extra_edge(PEB p);
	edge_descriptor compute_parent(edge_descriptor e);
};

#endif
