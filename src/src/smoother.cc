#include "smoother.h"

smoother::smoother(splice_graph &g)
	: gr(g)
{
	env = new GRBEnv();
	model = new GRBModel(*env);
}

smoother::~smoother()
{
	delete model;
	delete env;
}

int smoother::solve()
{
	build_edge_map();

	add_vertex_weight_variables();
	add_vertex_error_variables();
	add_edge_weight_variables();
	add_edge_error_variables();

	add_vertex_error_constraints();
	add_edge_error_constraints();

	add_conservation_constraints();

	model->getEnv().set(GRB_IntParam_OutputFlag, 0);

	model->optimize();

	update();

	return 0;
}

int smoother::build_edge_map()
{
	int index = 0;
	edge_iterator it1, it2;
	for(tie(it1, it2) = edges(gr); it1 != it2; it1++)
	{
		e2i.insert(PEI(*it1, index));
		i2e.push_back(*it1);
		index++;
	}
	return 0;
}

int smoother::add_vertex_weight_variables()
{
	vnwt.clear();
	for(int i = 0; i < num_vertices(gr); i++)
	{
		GRBVar var = model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
		vnwt.push_back(var);
	}
	model->update();
	return 0;
}

int smoother::add_vertex_error_variables()
{
	verr.clear();
	for(int i = 0; i < num_vertices(gr); i++)
	{
		double sd = get(get(vertex_stddev, gr), i);
		GRBVar var = model->addVar(0, GRB_INFINITY, 1.0 / sd, GRB_CONTINUOUS);
		verr.push_back(var);
	}
	model->update();
	return 0;
}

int smoother::add_edge_weight_variables()
{
	enwt.clear();
	for(int i = 0; i < num_edges(gr); i++)
	{
		GRBVar var = model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
		enwt.push_back(var);
	}
	model->update();
	return 0;
}

int smoother::add_edge_error_variables()
{
	eerr.clear();
	for(int i = 0; i < i2e.size(); i++)
	{
		double sd = get(get(edge_stddev, gr), i2e[i]);
		GRBVar var = model->addVar(0, GRB_INFINITY, 1.0 / sd, GRB_CONTINUOUS);
		eerr.push_back(var);
	}
	model->update();
	return 0;
}

int smoother::add_vertex_error_constraints()
{
	for(int i = 0; i < num_vertices(gr); i++)
	{
		double w = get(get(vertex_weight, gr), i);
		model->addConstr(vnwt[i], GRB_LESS_EQUAL, w + verr[i]);
		model->addConstr(vnwt[i], GRB_GREATER_EQUAL, w - verr[i]);
	}
	return 0;
}

int smoother::add_edge_error_constraints()
{
	for(int i = 0; i < num_edges(gr); i++)
	{
		double w = get(get(edge_weight, gr), i2e[i]);
		model->addConstr(enwt[i], GRB_LESS_EQUAL, w + eerr[i]);
		model->addConstr(enwt[i], GRB_GREATER_EQUAL, w - eerr[i]);
	}
	return 0;
}

int smoother::add_conservation_constraints()
{
	for(int i = 1; i < num_vertices(gr) - 1; i++)
	{
		int xi = in_degree(i, gr);
		int yi = out_degree(i, gr);
		if(xi == 0 && yi == 0) continue;

		assert(xi > 0 && yi > 0);

		in_edge_iterator i1, i2;
		GRBLinExpr ei;
		for(tie(i1, i2) = in_edges(i, gr); i1 != i2; i1++)
		{
			int e = e2i[*i1];
			assert(e >= 0 && e < num_edges(gr));
			ei += enwt[e];
		}

		out_edge_iterator o1, o2;
		GRBLinExpr eo;
		for(tie(o1, o2) = out_edges(i, gr); o1 != o2; o1++)
		{
			int e = e2i[*o1];
			assert(e >= 0 && e < num_edges(gr));
			eo += enwt[e];
		}

		model->addConstr(vnwt[i], GRB_EQUAL, ei);
		model->addConstr(vnwt[i], GRB_EQUAL, eo);
	}
	return 0;
}

int smoother::update()
{
	for(int i = 0; i < vnwt.size(); i++)
	{
		double w1 = vnwt[i].get(GRB_DoubleAttr_X);
		double w0 = get(get(vertex_weight, gr), i);
		put(get(vertex_weight, gr), i, w1);
	}

	for(int i = 0; i < enwt.size(); i++)
	{
		double w1 = enwt[i].get(GRB_DoubleAttr_X);
		double w0 = get(get(edge_weight, gr), i2e[i]);
		put(get(edge_weight, gr), i2e[i], w1);
	}

	return 0;
}
