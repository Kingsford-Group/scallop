#ifndef __LPSOLVER_H__
#define __LPSOLVER_H__

#include "gurobi_c++.h"
#include "dgraph.h"

class lpsolver
{
public:
	lpsolver(dgraph &g);
	virtual ~lpsolver();

private:
	dgraph &gr;						// splice graph
	vector<edge_descriptor> i2e;	// edge map
	MEI e2i;						// edge map

	vector<GRBVar> vnwt;			// new weights for nodes
	vector<GRBVar> enwt;			// new weights for 
	vector<GRBVar> verr;			// error for nodes
	vector<GRBVar> eerr;			// error for edges

	GRBModel * model;
	GRBEnv * env;

public:
	int solve();

private:
	int build_edge_map();
	int add_vertex_weight_variables();
	int add_vertex_error_variables();
	int add_edge_weight_variables();
	int add_edge_error_variables();

	int add_vertex_error_constraints();
	int add_edge_error_constraints();

	int add_conservation_constraints();

	int update();
};

#endif
