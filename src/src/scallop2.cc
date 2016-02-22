#include "scallop2.h"
#include <cstdio>

scallop2::scallop2(splice_graph &gr)
	: assembler(gr)
{
	GRBEnv * env = new GRBEnv();
	GRBModel * model = new GRBModel(*env);
}

scallop2::~scallop2()
{
	delete model;
	delete env;
}

int scallop2::assemble()
{
	smooth_weights();
	iterate();
	return 0;
}

int scallop2::iterate()
{
	init();
	while(true)
	{
		assign_weights();
		path p = compute_maximum_forward_path();
		cpaths.push_back(p);
		update_counters(p);
		update_lpsolver(p);
		double w = optimize();
		if(w <= cabd) break;
		cabd = w;
	}
	return 0;
}

int scallop2::init()
{
	// save the current weights
	get_edge_weights(gr, ewrt);
	assert(ewrt.size() == num_edges(gr));

	// init counters and current paths
	cabd = 0;
	cpaths.clear();
	ecnt.clear();
	edge_iterator it1, it2;
	for(tie(it1, it2) = edges(gr); it1 != it2; it1++)
	{
		ecnt.insert(PEI(*it1, 0));
	}

	// GRBLinExpr
	vars.clear();
	for(tie(it1, it2) = edges(gr); it1 != it2; it1++)
	{
		eexprs.insert(PEG(*it1, 0));
	}

	return 0;
}

int scallop2::assign_weights()
{
	edge_iterator it1, it2;
	for(tie(it1, it2) = edges(gr); it1 != it2; it1++)
	{
		assert(ecnt.find(*it1) != ecnt.end());
		assert(ewrt.find(*it1) != ewrt.end());
		double w = ewrt[*it1] / (1.0 + ecnt[*it1]);
		put(get(edge_weight, gr), *it1, w);
	}
	return 0;
}

int scallop2::update_counters(const path &p)
{
	for(int i = 0; i < p.v.size() - 1; i++)
	{
		PEB e = edge(p.v[i], p.v[i + 1], gr);
		assert(e.second == true);
		ecnt[e.first]++;
	}
	return 0;
}

int scallop2::update_lpsolver(const path &p)
{
	// add new variable for p
	GRBVar var = model->addVar(0, GRB_INFINITY, 1,  GRB_CONTINUOUS);
	vars.push_back(var);
	model->update();

	// update GRBLinExpr and constraints
	for(int i = 0; i < p.v.size() - 1; i++)
	{
		PEB e = edge(p.v[i], p.v[i + 1], gr);
		assert(e.second == true);
		
		eexprs[e.first] += var;
		model->addConstr(eexprs[e.first], GRB_LESS_EQUAL, ewrt[e.first]);
	}
	return 0;
}

double scallop2::optimize()
{
	model->optimize();
	for(int i = 0; i < vars.size(); i++)
	{
		double w = vars[i].get(GRB_DoubleAttr_X);
		cpaths[i].abd = w;
	}
	return model->get(GRB_DoubleAttr_ObjVal);
}

int scallop2::collect()
{
	// collect paths
	for(int i = 0; i < cpaths.size(); i++)
	{
		paths.push_back(cpaths[i]);
	}

	// update graph
	for(int i = 0; i < cpaths.size(); i++)
	{
		decrease_path(cpaths[i]);
	}
	return 0;
}
