/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __ROUTER_H__
#define __ROUTER_H__

#include <vector>
#include "util.h"
#include "splice_graph.h"
#include "equation.h"
#include "undirected_graph.h"
#include "hyper_set.h"

typedef pair<int, double> PID;
typedef map<int, double> MID;
typedef pair<PI, double> PPID;
typedef map<PI, double> MPID;

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
	MPID pe2w;					// decompose results (for pairs of edges)

#ifdef USECLP
	MID se2w;
#endif

public:
	int classify();												// compute status
	int build();												// give solution

	// init
	int build_indices();										// build u2e and e2u
	int build_bipartite_graph();								// build bipartite graph
	vector<double> compute_balanced_weights();					// balanced weights

	// decompose splitable vertex
	int split();												// for splitable vertices

	// decompose unsplitable vertex with greedy algorithm
	int thread();												// for unsplitable vertices
	int thread_isolate1(int k, vector<double> &vw);
	int thread_isolate2(int k, vector<double> &vw);
	bool thread_leaf(vector<double> &vw);
	bool thread_turn(vector<double> &vw);

#ifdef USECLP
	// decompose unsplitable vertex with LP 
	int lpsolve();
	int extend_bipartite_graph_max();							// extended graph
	int extend_bipartite_graph_all();							// extended graph
	int build_maximum_spanning_tree();							// make ug a (maximum) spanning tree
	int decompose0_clp();										// solve LP with CLP
	int decompose1_clp();										// solve LP with CLP
	int decompose2_clp();										// solve LP with CLP
#endif

	// print and stats
	int print();
	int stats();
};

bool compare_edge_weight(const PED &x, const PED &y);

#endif
