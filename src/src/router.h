#ifndef __ROUTER_H__
#define __ROUTER_H__

#include <vector>
#include "util.h"
#include "splice_graph.h"
#include "equation.h"
#include "undirected_graph.h"

using namespace std;

class router
{
public:
	router(int r, splice_graph &g, MEI &ei, VE &ie);
	router(int r, splice_graph &g, MEI &ei, VE &ie, const vector<PI> &p);
	router& operator=(const router &rt);

public:
	int root;					// central vertex
	splice_graph &gr;			// reference splice graph
	MEI &e2i;					// reference map of edge to index
	VE &i2e;					// reference map of index to edge
	vector<PI> routes;			// pairs of connections

	MI e2u;						// edge to index
	vector<int> u2e;			// index to edge
	undirected_graph ug;		// bipartite graph

	bool phasing;				// only use routes to divide
	double ratio;				// worst ratio
	vector<equation> eqns;		// divide results

public:
	// recompute everything
	int build(bool phasing);

	int build_indices();					// build u2e and e2u
	int build_bipartite_graph();			// build bipartite graph
	int divide();							// split *this, fill eqns
	int add_single_equation();				// cannot be divided
	int run_subsetsum();					// use subsetsum4
	int run_ilp1();							// with multiplier ratio
	int run_ilp2();							// with difference ratio

	// print and stats
	int print() const;
	int stats();
};

#endif
