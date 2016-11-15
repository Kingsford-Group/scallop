#ifndef __ROUTER_H__
#define __ROUTER_H__

#include <vector>
#include "util.h"
#include "splice_graph.h"
#include "equation.h"
#include "undirected_graph.h"

using namespace std;

typedef pair<PI, double> PPID;

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

	int status;					// splitable or not
	double ratio;				// worst ratio
	double obj1;				// increase of objective for splitable vertex
	double obj2;				// increase of objective for insplitable vertex
	vector<equation> eqns;		// decomposotion result for splitable vertex
	vector<PPID> vpi;			// decomposition result for insplitable vertex

public:
	// recompute everything
	int build();

	int build_indices();					// build u2e and e2u
	int build_bipartite_graph();			// build bipartite graph
	int classify();							// classify, splitable/insplitable
	int add_single_equation();				// cannot be divided
	int split();							// use subsetsum4

	// print and stats
	int print() const;
	int stats();
};

#endif
