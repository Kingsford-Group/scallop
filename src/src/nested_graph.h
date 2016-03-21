#ifndef __NESTED_GRAPH_H__
#define __NESTED_GRAPH_H__

#include "directed_graph.h"
#include "undirected_graph.h"

#include <map>
#include <cassert>

using namespace std;

typedef pair<int, int> PI;

class nested_graph : public directed_graph
{
public:
	nested_graph();
	virtual ~nested_graph();

public:
	vector<int> order;		// partial order of vertex
	vector<PI> partners;	// in/out partners for each vertex
	VE i2e;					// edge map
	MEI e2i;				// edge map
	vector<int>	parents;	// parent for each edge
	vector<bool> docking;	// whether this edge can be moved to boundary

public:
	int solve(directed_graph &gr);
	int draw(const string &file);
	int print();

private:
	int init(directed_graph &gr);
	int build_nests(directed_graph &gr);
	int build_partial_order();

	int build_partners(directed_graph &gr);
	int build_parents();
	int build_parents(int x);

	int build_docking();
	vector<int> bfs_docking_forward(int s, int p);
	vector<int> bfs_docking_backward(int t, int p);

	int compute_lca(int x, int y, vector<int> &xv, vector<int> &yv);
	bool dock(int e, int p);

	bool link(int xs, int xt, int ys, int yt, vector<PI> &p);
	PEB add_extra_edge(int s, int t);
	int remove_extra_edge(PEB p);
	int compute_parent(edge_descriptor e);
	int test_linking(directed_graph &gr);
};

#endif
