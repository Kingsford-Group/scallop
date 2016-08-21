#include "router.h"
#include "util.h"
#include "gurobi_c++.h"
#include <cstdio>
#include <algorithm>
#include <set>

router::router(int r, splice_graph &g, MEI &ei, VE &ie)
	:root(r), gr(g), e2i(ei), i2e(ie)
{
	build_indices();
	pratio = -1;
	pvalue = -1;
}

int router::divide()
{
	if(gr.in_degree(root) <= 1) return 0;
	if(gr.out_degree(root) <= 1) return 0;

	run_ilp1();
	evaluate_partition();
	print();
	printf("\n");

	run_ilp2();
	evaluate_partition();
	print();
	printf("\n");

	return 0;
}

int router::run_ilp1()
{
	GRBEnv *env = new GRBEnv();
	GRBModel *model = new GRBModel(*env);

	// ratio variables
	GRBVar xvar = model->addVar(1, GRB_INFINITY, 1, GRB_CONTINUOUS);
	model->update();

	// partition and weights variables
	vector<GRBVar> bvars;
	vector<GRBVar> pvars;
	vector<GRBVar> qvars;
	for(int i = 0; i < e2u.size(); i++)
	{
		GRBVar bvar = model->addVar(0, 1, 0, GRB_BINARY);
		GRBVar pvar = model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
		GRBVar qvar = model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
		bvars.push_back(bvar);
		pvars.push_back(pvar);
		qvars.push_back(qvar);
	}
	model->update();

	// constraints from routes
	for(int i = 0; i < routes.size(); i++)
	{
		int u1 = e2u[routes[i].first];
		int u2 = e2u[routes[i].second];
		assert(u1 >= 0 && u1 < bvars.size());
		assert(u2 >= 0 && u2 < bvars.size());
		model->addConstr(bvars[u1], GRB_EQUAL, bvars[u2]);
	}

	// guarantee a non-trivial partition
	GRBLinExpr exp1 = 0, exp2 = 0, exp3 = 0, exp4 = 0;
	for(int i = 0; i < gr.in_degree(root); i++)
	{
		exp1 += bvars[i];
		exp2 += (1 - bvars[i]);
	}
	for(int i = gr.in_degree(root); i < gr.degree(root); i++)
	{
		exp3 += bvars[i];
		exp4 += (1 - bvars[i]);
	}
	model->addConstr(exp1, GRB_GREATER_EQUAL, 1);
	model->addConstr(exp2, GRB_GREATER_EQUAL, 1);
	model->addConstr(exp3, GRB_GREATER_EQUAL, 1);
	model->addConstr(exp4, GRB_GREATER_EQUAL, 1);

	// lower bounds and upper bounds
	double ubound = 0.0;
	vector<double> weights;
	for(int i = 0; i < u2e.size(); i++)
	{
		edge_descriptor e = i2e[u2e[i]];
		double w = gr.get_edge_weight(e);
		weights.push_back(w);
		ubound += w;
	}
	for(int i = 0; i < bvars.size(); i++)
	{
		//printf("weight at %d = %.2lf, ubound = %.2lf\n", i, weights[i], ubound);
		model->addConstr(pvars[i], GRB_GREATER_EQUAL, bvars[i] * weights[i]);
		model->addConstr(pvars[i], GRB_LESS_EQUAL, bvars[i] * ubound);
		model->addConstr(pvars[i], GRB_LESS_EQUAL, xvar * weights[i]);

		model->addConstr(qvars[i], GRB_GREATER_EQUAL, (1 - bvars[i]) * weights[i]);
		model->addConstr(qvars[i], GRB_LESS_EQUAL, (1 - bvars[i]) * ubound);
		model->addConstr(qvars[i], GRB_LESS_EQUAL, xvar * weights[i]);
	}

	// balance
	exp1 = exp2 = exp3 = exp4 = 0;
	for(int i = 0; i < gr.in_degree(root); i++)
	{
		exp1 += pvars[i];
		exp2 += qvars[i];
	}
	for(int i = gr.in_degree(root); i < gr.degree(root); i++)
	{
		exp3 += pvars[i];
		exp4 += qvars[i];
	}
	model->addConstr(exp1, GRB_EQUAL, exp3);
	model->addConstr(exp2, GRB_EQUAL, exp4);

	// set objective and optimize
	GRBLinExpr obj = xvar;
	model->setObjective(obj, GRB_MINIMIZE);
	model->getEnv().set(GRB_IntParam_OutputFlag, 0);
	model->update();

	model->optimize();

	int f = model->get(GRB_IntAttr_Status);

	pv1.clear();
	pv2.clear();
	qv1.clear();
	qv2.clear();
	if(f == GRB_OPTIMAL)
	{
		for(int i = 0; i < bvars.size(); i++)
		{
			int b = (int)(bvars[i].get(GRB_DoubleAttr_X));
			double p = pvars[i].get(GRB_DoubleAttr_X);
			double q = qvars[i].get(GRB_DoubleAttr_X);

			if(b == 1) assert(q <= 0.001);
			if(b == 0) assert(p <= 0.001);
			if(b == 1 && i < gr.in_degree(root)) pv1.push_back(u2e[i]);
			if(b == 1 && i >=gr.in_degree(root)) pv2.push_back(u2e[i]);
			if(b == 0 && i < gr.in_degree(root)) qv1.push_back(u2e[i]);
			if(b == 0 && i >=gr.in_degree(root)) qv2.push_back(u2e[i]);
		}
	}

	delete model;
	delete env;
}

int router::run_ilp2()
{
	GRBEnv *env = new GRBEnv();
	GRBModel *model = new GRBModel(*env);

	// ratio variables
	GRBVar xvar = model->addVar(1, GRB_INFINITY, 1, GRB_CONTINUOUS);
	model->update();

	// partition and weights variables
	vector<GRBVar> bvars;
	for(int i = 0; i < e2u.size(); i++)
	{
		GRBVar bvar = model->addVar(0, 1, 0, GRB_BINARY);
		bvars.push_back(bvar);
	}
	model->update();

	// constraints from routes
	for(int i = 0; i < routes.size(); i++)
	{
		int u1 = e2u[routes[i].first];
		int u2 = e2u[routes[i].second];
		assert(u1 >= 0 && u1 < bvars.size());
		assert(u2 >= 0 && u2 < bvars.size());
		model->addConstr(bvars[u1], GRB_EQUAL, bvars[u2]);
	}

	// guarantee a non-trivial partition
	GRBLinExpr exp1 = 0, exp2 = 0, exp3 = 0, exp4 = 0;
	for(int i = 0; i < gr.in_degree(root); i++)
	{
		exp1 += bvars[i];
		exp2 += (1 - bvars[i]);
	}
	for(int i = gr.in_degree(root); i < gr.degree(root); i++)
	{
		exp3 += bvars[i];
		exp4 += (1 - bvars[i]);
	}
	model->addConstr(exp1, GRB_GREATER_EQUAL, 1);
	model->addConstr(exp2, GRB_GREATER_EQUAL, 1);
	model->addConstr(exp3, GRB_GREATER_EQUAL, 1);
	model->addConstr(exp4, GRB_GREATER_EQUAL, 1);

	// lower bounds and upper bounds
	vector<double> weights;
	for(int i = 0; i < u2e.size(); i++)
	{
		edge_descriptor e = i2e[u2e[i]];
		double w = gr.get_edge_weight(e);
		weights.push_back(w);
	}

	// balance
	GRBLinExpr pexp1 = 0, pexp2 = 0, qexp1 = 0, qexp2 = 0;
	for(int i = 0; i < gr.in_degree(root); i++)
	{
		pexp1 += bvars[i] * weights[i];
		qexp1 += (1 - bvars[i]) * weights[i];
	}
	for(int i = gr.in_degree(root); i < gr.degree(root); i++)
	{
		pexp2 += bvars[i] * weights[i];
		qexp2 += (1 - bvars[i]) * weights[i];
	}
	model->addConstr(pexp1, GRB_GREATER_EQUAL, pexp2 - xvar);
	model->addConstr(pexp1, GRB_LESS_EQUAL, pexp2 + xvar);
	model->addConstr(qexp1, GRB_GREATER_EQUAL, qexp2 - xvar);
	model->addConstr(qexp1, GRB_LESS_EQUAL, qexp2 + xvar);

	// set objective and optimize
	GRBLinExpr obj = xvar;
	model->setObjective(obj, GRB_MINIMIZE);
	model->getEnv().set(GRB_IntParam_OutputFlag, 0);
	model->update();

	model->optimize();

	int f = model->get(GRB_IntAttr_Status);

	pv1.clear();
	pv2.clear();
	qv1.clear();
	qv2.clear();
	if(f == GRB_OPTIMAL)
	{
		for(int i = 0; i < bvars.size(); i++)
		{
			int b = (int)(bvars[i].get(GRB_DoubleAttr_X));
			if(b == 1 && i < gr.in_degree(root)) pv1.push_back(u2e[i]);
			if(b == 1 && i >=gr.in_degree(root)) pv2.push_back(u2e[i]);
			if(b == 0 && i < gr.in_degree(root)) qv1.push_back(u2e[i]);
			if(b == 0 && i >=gr.in_degree(root)) qv2.push_back(u2e[i]);
		}
	}

	delete model;
	delete env;
}

int router::evaluate_partition()
{
	if(pv1.size() <= 0 || pv2.size() <= 0) return 0;
	if(qv1.size() <= 0 || qv2.size() <= 0) return 0;

	double pw1 = 0, pw2 = 0;
	double qw1 = 0, qw2 = 0;

	for(int i = 0; i < pv1.size(); i++) pw1 += gr.get_edge_weight(i2e[pv1[i]]);
	for(int i = 0; i < pv2.size(); i++) pw2 += gr.get_edge_weight(i2e[pv2[i]]);
	for(int i = 0; i < qv1.size(); i++) qw1 += gr.get_edge_weight(i2e[qv1[i]]);
	for(int i = 0; i < qv2.size(); i++) qw2 += gr.get_edge_weight(i2e[qv2[i]]);

	double pr = (pw1 > pw2) ? (pw1 / pw2) : (pw2 / pw1);
	double qr = (qw1 > qw2) ? (qw1 / qw2) : (qw2 / qw1);
	pratio = (pr < qr) ? qr : pr;

	double pw = fabs(pw1 - pw2);
	double qw = fabs(qw1 - qw2);
	pvalue = (pw < qw) ? qw : pw;

	return 0;
}

int router::build_indices()
{
	e2u.clear();
	u2e.clear();

	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.in_edges(root); it1 != it2; it1++)
	{
		int e = e2i[*it1];
		e2u.insert(PI(e, e2u.size()));
		u2e.push_back(e);
	}
	for(tie(it1, it2) = gr.out_edges(root); it1 != it2; it1++)
	{
		int e = e2i[*it1];
		e2u.insert(PI(e, e2u.size()));
		u2e.push_back(e);
	}

	return 0;
}

int router::add_route(const PI &p, double c)
{
	for(int i = 0; i < routes.size(); i++)
	{
		if(routes[i] == p)
		{
			counts[i] += c;
			return 0;
		}
	}

	routes.push_back(p);
	counts.push_back(c);
	return 0;
}

int router::remove_route(const PI &p)
{
	for(int i = 0; i < routes.size(); i++)
	{
		if(routes[i] == p)
		{
			routes.erase(routes.begin() + i);
			counts.erase(counts.begin() + i);
			return 0;
		}
	}
	return 0;
}

int router::replace_in_edge(int ex, int ey)
{
	for(int i = 0; i < routes.size(); i++)
	{
		if(routes[i].first == ex) routes[i].first = ey;
	}
	return 0;
}

int router::replace_out_edge(int ex, int ey)
{
	for(int i = 0; i < routes.size(); i++)
	{
		if(routes[i].second == ex) routes[i].second = ey;
	}
	return 0;
}

int router::split_in_edge(int ex, int ey, double r)
{
	assert(r >= 0 && r <= 1.0);
	if(ex == ey) return 0;
	int n = routes.size();
	for(int i = 0; i < n; i++)
	{
		if(routes[i].first == ex)
		{
			add_route(PI(ey, routes[i].second), (1.0 - r) * counts[i]);
			counts[i] *= r;
		}
	}
	return 0;
}

int router::split_out_edge(int ex, int ey, double r)
{
	assert(r >= 0 && r <= 1.0);
	if(ex == ey) return 0;
	int n = routes.size();
	for(int i = 0; i < n; i++)
	{
		if(routes[i].second == ex)
		{
			add_route(PI(routes[i].first, ey), (1.0 - r) * counts[i]);
			counts[i] *= r;
		}
	}
	return 0;
}

int router::remove_in_edges(const vector<int> &v)
{
	set<int> s(v.begin(), v.end());

	vector<PI> vv;
	vector<double> cc;
	for(int i = 0; i < routes.size(); i++)
	{
		int x = routes[i].first;
		if(s.find(x) != s.end()) continue;
		vv.push_back(routes[i]);
		cc.push_back(counts[i]);
	}

	routes = vv;
	counts = cc;

	return 0;
}

int router::remove_out_edges(const vector<int> &v)
{
	set<int> s(v.begin(), v.end());

	vector<PI> vv;
	vector<double> cc;
	for(int i = 0; i < routes.size(); i++)
	{
		int x = routes[i].second;
		if(s.find(x) != s.end()) continue;
		vv.push_back(routes[i]);
		cc.push_back(counts[i]);
	}

	routes = vv;
	counts = cc;

	return 0;
}

int router::remove_in_edge(int x)
{
	vector<int> v;
	v.push_back(x);
	return remove_in_edges(v);
}

int router::remove_out_edge(int x)
{
	vector<int> v;
	v.push_back(x);
	return remove_out_edges(v);
}

double router::total_counts() const
{
	double s = 0;
	for(int i = 0; i < counts.size(); i++)
	{
		s += counts[i];
	}
	return s;
}

int router::print() const
{
	printf("router %d, #routes = %lu, total-counts = %.1lf\n", root, routes.size(), total_counts());
	for(int i = 0; i < routes.size(); i++)
	{
		printf(" route %d (%d, %d), count = %.1lf\n", i, routes[i].first, routes[i].second, counts[i]);
	}
	printf("divide vertex %d, pratio = %.2lf, pvalue = %.2lf, pv1 = ( ", root, pratio, pvalue);
	printv(pv1);
	printf("), pv2 = ( ");
	printv(pv2);
	printf("), qv1 = ( ");
	printv(qv1);
	printf("), qv2 = ( ");
	printv(qv2);
	printf(")\n");
	return 0;
}

