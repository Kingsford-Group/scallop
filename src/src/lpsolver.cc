#include "lpsolver.h"

lpsolver::lpsolver(const bbase &bb)
	: sgraph(bb)
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
