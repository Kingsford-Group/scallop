#include "lpsolver.h"

lpsolver::lpsolver(const dgraph &g)
	: gr(g)
{
	env = new GRBEnv();
	model = new GRBModel(*env);
}

lpsolver::~lpsolver()
{
	delete model;
	delete env;
}

int lpsolver::solve()
{
	return 0;
}
