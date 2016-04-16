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
	add_edge_weight_constraints();
	add_conservation_constraints();
	add_additional_constraints();

	model->getEnv().set(GRB_IntParam_OutputFlag, 0);

	model->optimize();

	update();

	return 0;
}

int smoother::add_equation(const VE &x, const VE &y)
{
	xeq.push_back(x);
	yeq.push_back(y);
	return 0;
}

int smoother::build_edge_map()
{
	int index = 0;
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
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
	for(int i = 0; i < gr.num_vertices(); i++)
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
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		double sd = gr.get_vertex_stddev(i);
		GRBVar var = model->addVar(0, GRB_INFINITY, 1.0 / sd, GRB_CONTINUOUS);
		verr.push_back(var);
	}
	model->update();
	return 0;
}

int smoother::add_edge_weight_variables()
{
	enwt.clear();
	for(int i = 0; i < gr.num_edges(); i++)
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
		double sd = gr.get_edge_stddev(i2e[i]);
		GRBVar var = model->addVar(0, GRB_INFINITY, 1.0 / sd, GRB_CONTINUOUS);
		eerr.push_back(var);
	}
	model->update();
	return 0;
}

int smoother::add_vertex_error_constraints()
{
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		double w = gr.get_vertex_weight(i);
		model->addConstr(vnwt[i], GRB_LESS_EQUAL, w + verr[i]);
		model->addConstr(vnwt[i], GRB_GREATER_EQUAL, w - verr[i]);
	}
	return 0;
}

int smoother::add_edge_error_constraints()
{
	for(int i = 0; i < gr.num_edges(); i++)
	{
		double w = gr.get_edge_weight(i2e[i]);
		model->addConstr(enwt[i], GRB_LESS_EQUAL, w + eerr[i]);
		model->addConstr(enwt[i], GRB_GREATER_EQUAL, w - eerr[i]);
	}
	return 0;
}

int smoother::add_edge_weight_constraints()
{
	for(int i = 0; i < gr.num_edges(); i++)
	{
		model->addConstr(enwt[i], GRB_GREATER_EQUAL, 1.0);
	}
	return 0;
}

int smoother::add_conservation_constraints()
{
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		int xi = gr.in_degree(i);
		int yi = gr.out_degree(i);
		if(xi == 0 && yi == 0) continue;

		assert(xi > 0 && yi > 0);

		edge_iterator i1, i2;
		GRBLinExpr ei;
		for(tie(i1, i2) = gr.in_edges(i); i1 != i2; i1++)
		{
			int e = e2i[*i1];
			assert(e >= 0 && e < gr.num_edges());
			ei += enwt[e];
		}

		edge_iterator o1, o2;
		GRBLinExpr eo;
		for(tie(o1, o2) = gr.out_edges(i); o1 != o2; o1++)
		{
			int e = e2i[*o1];
			assert(e >= 0 && e < gr.num_edges());
			eo += enwt[e];
		}

		model->addConstr(vnwt[i], GRB_EQUAL, ei);
		model->addConstr(vnwt[i], GRB_EQUAL, eo);
	}
	return 0;
}

int smoother::add_additional_constraints()
{
	assert(xeq.size() == yeq.size());
	for(int i = 0; i < xeq.size(); i++)
	{
		GRBLinExpr ei;
		for(int k = 0; k < xeq[i].size(); k++)
		{
			int e = e2i[xeq[i][k]];
			assert(e >= 0 && e < gr.num_edges());
			ei += enwt[e];
		}

		GRBLinExpr eo;
		for(int k = 0; k < yeq[i].size(); k++)
		{
			int e = e2i[yeq[i][k]];
			assert(e >= 0 && e < gr.num_edges());
			eo += enwt[e];
		}

		model->addConstr(ei, GRB_EQUAL, eo);
	}
	return 0;
}

int smoother::update()
{
	for(int i = 0; i < vnwt.size(); i++)
	{
		double w1 = vnwt[i].get(GRB_DoubleAttr_X);
		double w0 = gr.get_vertex_weight(i);
		gr.set_vertex_weight(i, w1);
	}

	for(int i = 0; i < enwt.size(); i++)
	{
		double w1 = enwt[i].get(GRB_DoubleAttr_X);
		double w0 = gr.get_edge_weight(i2e[i]);
		gr.set_edge_weight(i2e[i], w1);
	}

	return 0;
}
