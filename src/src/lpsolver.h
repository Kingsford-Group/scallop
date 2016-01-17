#ifndef __LPSOLVER_H__
#define __LPSOLVER_H__

#include "gurobi_c++.h"
#include "sgraph.h"

class lpsolver : public sgraph
{
public:
	lpsolver(const bbase &bb);
	virtual ~lpsolver();

public:
	// members for ILP
	vector<GRBVar> vars;					// gene variables
	GRBModel * model;
	GRBEnv * env;

public:
	int solve();
};

#endif
