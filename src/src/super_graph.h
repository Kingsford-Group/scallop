#ifndef __SUPER_GRAPH_H__
#define __SUPER_GRAPH_H__

#include "hyper_set.h"
#include "undirected_graph.h"
#include "splice_graph.h"
#include "util.h"

#include <map>
#include <cassert>

using namespace std;

class super_graph
{
public:
	super_graph(const splice_graph &gr, const hyper_set &hs);
	virtual ~super_graph();

public:
	splice_graph root;			// splice graph
	hyper_set hyper;			// hyper set
	vector<splice_graph> subs;	// sub-graphs
	vector<hyper_set> hss;		// sub-hyper-set

private:
	undirected_graph ug;		// graph without edges to s and t
	map<int, PI> a2b;			// vertex map from gr to subgraphs
	map<PI, int> b2a;			// vertex map from subgraphs to gr

public:
	int build();
	int print();
	int get_root_vertex(int sub, int x) const;
	vector<int> get_root_vertices(int sub, const vector<int> &x) const;

private:
	int build_undirected_graph();
	int split_splice_graph();
	int build_single_splice_graph(splice_graph &gr, hyper_set &hs, const set<int> &v, int index);

	// analysis the structure
	bool remove_single_read(splice_graph &gr);
	int build_cuts(splice_graph &gr);
	int build_maximum_path_graph(splice_graph &gr, undirected_graph &mg);
	double compute_maximum_path1(splice_graph &gr, int s, int &t);
	double compute_maximum_path2(splice_graph &gr, int t, int &s);
};

#endif
