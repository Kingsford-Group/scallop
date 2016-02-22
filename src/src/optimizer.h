#ifndef __LPSOLVER_H__
#define __LPSOLVER_H__

#include "gurobi_c++.h"
#include "splice_graph.h"
#include "path.h"

class optimizer
{
public:
	optimizer(const splice_graph &g, vector<path> &paths);
	virtual ~optimizer();

private:
	const splice_graph &gr;			// splice graph
	vector<path> &paths;			// paths to be optimized
	vector<edge_descriptor> i2e;	// edge map
	MEI e2i;						// edge map

	vector<GRBVar> pnwt;			// new weights for paths

	GRBModel * model;
	GRBEnv * env;

public:
	int solve();

private:
	int build_edge_map();
	int add_path_weight_variables();

	int add_vertex_weight_constraints();
	int add_edge_weight_constraints();

	int update();
};

#endif
