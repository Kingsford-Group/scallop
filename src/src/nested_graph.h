#ifndef __NESTED_GRAPH_H__
#define __NESTED_GRAPH_H__

#include "directed_graph.h"
#include "undirected_graph.h"

#include <map>
#include <cassert>
#include <string>

using namespace std;

typedef pair<int, int> PI;

class nested_graph : public directed_graph
{
public:
	nested_graph();
	virtual ~nested_graph();

public:
	VE i2e;					// edge map
	MEI e2i;				// edge map
	vector<int> order;		// partial order of vertex
	vector<PI> partners;	// in/out partners for each vertex
	vector<int>	parents;	// parent for each edge

public:
	int clear();
	int build(directed_graph &gr);
	int draw(const string &file);
	bool link(int xs, int xt, int ys, int yt, vector<PI> &xp, vector<PI> &yp);

private:
	int init(directed_graph &gr);
	int build_nests(directed_graph &gr);
	int build_partial_order();
	int build_partners(directed_graph &gr);
	int build_parents();
	int build_parents(int x);

	bool bfs_docking_forward(int e, int t, int p, vector<int> &v);
	bool bfs_docking_backward(int e, int s, int p, vector<int> &v);

	int compute_lca(int x, int y, vector<int> &xv, vector<int> &yv);
	int dock(int e, int p, vector<PI> &v);

	PEB add_extra_edge(int s, int t);
	int remove_extra_edge(PEB p);
	int compute_parent(edge_descriptor e);
	
	int test_linking(directed_graph &gr);
	int print();
};

#endif
