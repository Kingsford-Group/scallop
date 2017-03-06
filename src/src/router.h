#ifndef __ROUTER_H__
#define __ROUTER_H__

#include <vector>
#include "util.h"
#include "splice_graph.h"
#include "equation.h"
#include "undirected_graph.h"
#include "hyper_set.h"

typedef pair<PI, double> PPID;

using namespace std;

class router
{
public:
	router(int r, splice_graph &g, MEI &ei, VE &ie);
	router(int r, splice_graph &g, MEI &ei, VE &ie, const MPII &mpi);
	router& operator=(const router &rt);

public:
	int root;					// central vertex
	splice_graph &gr;			// reference splice graph
	MEI &e2i;					// reference map of edge to index
	VE &i2e;					// reference map of index to edge
	vector<PI> routes;			// pairs of connections
	vector<int> counts;			// counts for routes

	MI e2u;						// edge to index
	vector<int> u2e;			// index to edge
	MED u2w;					// weights of edges of ug
	undirected_graph ug;		// bipartite graph

	int type;					// trivial, splitable, single, or multiple 
	int degree;					// level
	double ratio;				// worst ratio
	vector<equation> eqns;		// split results
	vector<PPID> vpi;			// decompose resutls

public:
	int classify();							// compute status
	int build();							// give solution

	int build_indices();					// build u2e and e2u
	int build_bipartite_graph();			// build bipartite graph
	int build_maximum_spanning_tree();	// make ug a (maximum) spanning tree
	int complete();							// complete graph
	int split();							// split
	int decompose1();						// decompose with quadratic linear programming
	int decompose2();						// decompose with linear system
	PI filter_hyper_edge();					// try to filter hyper-edge
	PI filter_small_hyper_edge();			// hyper-edge w.r.t. the smallest edge
	PI filter_cycle_hyper_edge();			// hyper-edge w.r.t. any cycle

	// print and stats
	int print() const;
	int stats();
};

#endif
