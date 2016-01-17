#ifndef __SGRAPH_H__
#define __SGRAPH_H__

#include "bundle.h"
#include "dgraph.h"

// class for splice graph
class sgraph: public bundle
{
public:
	sgraph(const bbase &bb);
	virtual ~sgraph();

public:
	dgraph gr;						// splice graph
	MEI e2b;						// edge-descriptor to bridge index
	MED e2w;						// edge-descriptor to weight

	vector<int> iv;					// pointers to each region (node)
	vector<int> ov;					// pointers from each region (node)
	vector<edge_descriptor> basis;	// representatives for the basis

public:
	int solve();
	int build_graph();
	int check();
	int print(int index);
	int draw(const string &file);

	int build_in_out_edges();
	int build_basis();

	PEB get_max_in_edge(int x);
	PEB get_max_out_edge(int x);
};

#endif
