#ifndef __SUPER_GRAPH_H__
#define __SUPER_GRAPH_H__

#include "hyper_edge.h"
#include "undirected_graph.h"
#include "splice_graph.h"
#include "util.h"

#include <map>
#include <cassert>

using namespace std;

class super_graph
{
public:
	super_graph(const splice_graph &gr);
	virtual ~super_graph();

public:
	splice_graph root;			// splice graph
	vector<splice_graph> subs;	// sub-graphs

private:
	undirected_graph mg;		// maximum path graph

	undirected_graph ug;		// graph without edges to s and t
	map<int, PI> a2b;			// vertex map from gr to subgraphs
	map<PI, int> b2a;			// vertex map from subgraphs to gr

public:
	int build();
	int print();
	int get_root_vertex(int sub, int x) const;
	vector<int> get_root_vertices(int sub, const vector<int> &x) const;

private:
	int build_maximum_path_graph();
	int build_undirected_graph();
	int split_splice_graph();
	int build_single_splice_graph(splice_graph &gr, const set<int> &v, int index);
};

#endif
