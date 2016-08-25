#ifndef __ROUTER_H__
#define __ROUTER_H__

#include <vector>
#include "util.h"
#include "splice_graph.h"
#include "equation.h"
#include "gurobi_c++.h"

using namespace std;

class router
{
public:
	router(int r, splice_graph &g, MEI &ei, VE &ie, GRBEnv *env);
	router(int r, splice_graph &g, MEI &ei, VE &ie);
	//router(const router &rt);
	router& operator=(const router &rt);

public:
	GRBEnv *env;				// GRB environment

	bool touched;				// whether it has been touched
	int root;					// central vertex
	splice_graph &gr;			// reference splice graph
	MEI &e2i;					// reference map of edge to index
	VE &i2e;					// reference map of index to edge

	vector<PI> routes;			// pairs of connections
	vector<double> counts;		// reads spanning each route
	MI e2u;						// edge to index
	vector<int> u2e;			// index to edge

	double ratio;				// worst ratio
	vector<equation> eqns;		// divide results

public:
	// recompute everything
	int update();

	int update_routes();			// remove false routes
	int build_indices();			// build u2e and e2u
	int divide();					// split *this, fill eqns
	int add_single_equation();		// cannot be divided
	int run_subsetsum();			// use subsetsum4
	int run_ilp1();					// with multiplier ratio
	int run_ilp2();					// with difference ratio
	int evaluate();					// evaluate eqns

	// modify routes
	int add_route(const PI &p, double c);
	int remove_route(const PI &p);
	int replace_in_edge(int ex, int ey);
	int replace_out_edge(int ex, int ey);
	int split_in_edge(int ex, int ey, double r);
	int split_out_edge(int ex, int ey, double r);
	int remove_in_edges(const vector<int> &v);
	int remove_out_edges(const vector<int> &v);
	int remove_in_edge(int x);
	int remove_out_edge(int x);

	// print and stats
	int print() const;
	double total_counts() const;
};

#endif
