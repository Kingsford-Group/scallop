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

	vector<double> bw;			// balanced weight
	double bratio;				// balanced ratio

	int status;					// splitable or not
	double ratio;				// split/decompose ratio
	double delta;				// decrease of the objective function
	vector<PPID> vpi;			// decomposition result for insplitable vertex
	vector<equation> eqns;		// decomposotion result for splitable vertex

public:
	// recompute everything
	int build();
	int solve();

	int print() const;
	int stats();

private:
	int build_indices();								// build u2e and e2u
	int build_bipartite_graph();						// build bipartite graph
	int build_balanced_weights();						// compute balanced weights
	int build_balanced_weights(const set<int> &fb);		// compute balanced weights
	int classify();										// classify, splitable/insplitable
	int decompose();									// for insplitable root
	int split();										// for splitable root
};

#endif
