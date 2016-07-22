#include "smoother.h"
#include "config.h"

#include <cmath>

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

int smoother::smooth_vertex(int i, double f1, double f2)
{
	try
	{
		build_edge_map(i);
		build_vertex_map(i);
		build_factor(i, f1, f2);

		add_vertex_weight_variables();
		add_vertex_error_variables();
		add_edge_weight_variables();
		add_edge_error_variables();

		add_vertex_error_constraints();
		add_edge_error_constraints();
		add_edge_weight_constraints();
		add_conservation_constraints();
		add_additional_constraints();

		set_objective();

		model->getEnv().set(GRB_IntParam_OutputFlag, 0);

		model->optimize();

		int f = model->get(GRB_IntAttr_Status);

		if (f != GRB_OPTIMAL) return -1;

		update();
	} 
	catch(GRBException e) 
	{
		cout << "GRB exception code = " << e.getErrorCode() << ", message = " << e.getMessage() << endl;
		return -1;
	} 
	catch(...) 
	{
		cout << "GRB exception" << endl;
		return -1;
	}

	return 0;
}


int smoother::smooth()
{
	try
	{
		build_edge_map();
		build_vertex_map();
		build_factors();

		add_vertex_weight_variables();
		add_vertex_error_variables();
		add_edge_weight_variables();
		add_edge_error_variables();

		add_vertex_error_constraints();
		add_edge_error_constraints();
		add_edge_weight_constraints();
		add_conservation_constraints();
		add_additional_constraints();

		set_objective();

		model->getEnv().set(GRB_IntParam_OutputFlag, 0);

		model->optimize();

		int f = model->get(GRB_IntAttr_Status);

		if (f != GRB_OPTIMAL) return -1;

		update();
	} 
	catch(GRBException e) 
	{
		cout << "GRB exception code = " << e.getErrorCode() << ", message = " << e.getMessage() << endl;
		return -1;
	} 
	catch(...) 
	{
		cout << "GRB exception" << endl;
		return -1;
	}

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
	i2e.clear();
	e2i.clear();
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		e2i.insert(PEI(*it1, index));
		i2e.push_back(*it1);
		index++;
	}
	return 0;
}

int smoother::build_edge_map(int i)
{
	i2e.clear();
	e2i.clear();
	int index = 0;
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.in_edges(i); it1 != it2; it1++)
	{
		e2i.insert(PEI(*it1, index));
		i2e.push_back(*it1);
		index++;
	}
	for(tie(it1, it2) = gr.out_edges(i); it1 != it2; it1++)
	{
		e2i.insert(PEI(*it1, index));
		i2e.push_back(*it1);
		index++;
	}
	return 0;
}

int smoother::build_vertex_map()
{
	v2i.clear();
	i2v.clear();
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		v2i.push_back(i);
		i2v.push_back(i);
	}
	return 0;
}

int smoother::build_vertex_map(int i)
{
	v2i.assign(gr.num_vertices(), -1);
	v2i[i] = 0;
	i2v.assign(1, 0);
	return 0;
}

int smoother::build_factors()
{
	sf.assign(i2v.size(), 1.0);
	tf.assign(i2v.size(), 1.0);
	return 0;
}

int smoother::build_factor(int i, double f1, double f2)
{
	sf.assign(1, f1);
	tf.assign(1, f2);
	return 0;
}


int smoother::add_vertex_weight_variables()
{
	vnwt.clear();
	for(int i = 0; i < i2v.size(); i++)
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
	for(int i = 0; i < i2v.size(); i++)
	{
		GRBVar var = model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
		verr.push_back(var);
	}
	model->update();
	return 0;
}

int smoother::add_edge_weight_variables()
{
	enwt.clear();
	for(int i = 0; i < i2e.size(); i++)
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
		GRBVar var = model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
		eerr.push_back(var);
	}
	model->update();
	return 0;
}

int smoother::add_vertex_error_constraints()
{
	for(int k = 0; k < i2v.size(); k++)
	{
		int i = i2v[k];
		if(gr.degree(i) == 0) continue;
		if(i == 0) continue;
		if(i == gr.num_vertices() - 1) continue;
		double w = gr.get_vertex_weight(i);
		model->addConstr(vnwt[i], GRB_LESS_EQUAL, w + verr[i]);
		model->addConstr(vnwt[i], GRB_GREATER_EQUAL, w - verr[i]);
	}
	return 0;
}

int smoother::add_edge_error_constraints()
{
	for(int i = 0; i < i2e.size(); i++)
	{
		double w = gr.get_edge_weight(i2e[i]);
		model->addConstr(enwt[i], GRB_LESS_EQUAL, w + eerr[i]);
		model->addConstr(enwt[i], GRB_GREATER_EQUAL, w - eerr[i]);
	}
	return 0;
}

int smoother::add_edge_weight_constraints()
{
	for(int i = 0; i < i2e.size(); i++)
	{
		model->addConstr(enwt[i], GRB_GREATER_EQUAL, 1.0);
	}
	return 0;
}

int smoother::add_conservation_constraints()
{
	for(int k = 0; k < i2v.size(); k++)
	{
		int i = i2v[k];
		if(i == 0) continue;
		if(i == gr.num_vertices() - 1) continue;

		//int xi = gr.in_degree(i);
		//int yi = gr.out_degree(i);
		//if(xi == 0 && yi == 0) continue;
		//assert(xi > 0 && yi > 0);

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

		model->addConstr(ei * sf[k], GRB_EQUAL, eo * tf[k]);
		model->addConstr(vnwt[k], GRB_EQUAL, (ei + eo) / 2.0);
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

int smoother::set_objective()
{
	GRBQuadExpr expr;
	for(int k = 0; k < i2v.size(); k++)
	{
		int i = i2v[k];
		if(i == 0) continue;
		if(i == gr.num_vertices() - 1) continue;
		double ll = gr.get_vertex_info(i).length * 1.0 / average_read_length;
		double wt = 1.0 + ll;
		expr += wt * verr[k] * verr[k];
	}

	for(int i = 0; i < i2e.size(); i++)
	{
		double ll = gr.get_edge_info(i2e[i]).length * 1.0 / average_read_length;
		double wt = 1.0 + ll;
		int s = i2e[i]->source();
		int t = i2e[i]->target();
		if(s == 0) wt = 0;
		if(t == gr.num_vertices() - 1) wt = 0;
		expr += wt * eerr[i] * eerr[i];
	}

	model->setObjective(expr, GRB_MINIMIZE);
	return 0;
}

int smoother::update()
{
	for(int k = 0; k < vnwt.size(); k++)
	{
		int i = i2v[k];
		double w1 = vnwt[k].get(GRB_DoubleAttr_X);
		gr.set_vertex_weight(i, w1);
	}

	for(int i = 0; i < enwt.size(); i++)
	{
		double w1 = enwt[i].get(GRB_DoubleAttr_X);
		edge_descriptor e = i2e[i];
		gr.set_edge_weight(i2e[i], w1);
	}
	return 0;
}
