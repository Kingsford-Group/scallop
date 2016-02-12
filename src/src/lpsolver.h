#ifndef __LPSOLVER_H__
#define __LPSOLVER_H__

#include "gurobi_c++.h"
#include "dgraph.h"

class lpsolver
{
public:
	lpsolver(const dgraph &g);
	virtual ~lpsolver();

public:
	const dgraph &gr;			// splice graph

	vector<GRBVar> vnwt;		// new weights for nodes
	vector<GRBVar> enwt;		// new weights for 
	vector<GRBVar> verr;		// error for nodes
	vector<GRBVar> eerr;		// error for edges

	GRBModel * model;
	GRBEnv * env;

public:
	int solve();
};

#endif
