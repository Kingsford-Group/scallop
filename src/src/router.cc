#include "router.h"
#include "util.h"
#include "gurobi_c++.h"
#include "smoother.h"
#include "subsetsum.h"
#include "config.h"
#include <cstdio>
#include <algorithm>
#include <set>

router::router(int r, splice_graph &g, MEI &ei, VE &ie)
	:root(r), gr(g), e2i(ei), i2e(ie), status(0)
{
}

router::router(int r, splice_graph &g, MEI &ei, VE &ie, const vector<PI> &p)
	:root(r), gr(g), e2i(ei), i2e(ie), routes(p), status(0)
{
}

router& router::operator=(const router &rt)
{
	root = rt.root;
	gr = rt.gr;
	e2i = rt.e2i;
	i2e = rt.i2e;

	routes = rt.routes;
	
	e2u = rt.e2u;
	u2e = rt.u2e;

	beta1 = rt.beta1;
	beta2 = rt.beta2;
	delta = rt.delta;
	status = rt.status;
	eqns = rt.eqns;
	vpi = rt.vpi;

	return (*this);
}

int router::build()
{
	assert(gr.in_degree(root) >= 1);
	assert(gr.out_degree(root) >= 1);

	delta = -1;
	status = -1;
	eqns.clear();
	vpi.clear();

	assert_balance();
	build_indices();
	build_bipartite_graph();
	compute_params();
	classify();

	if(status == SPLITABLE) split();
	else if(status == INSPLITABLE) decompose();

	return 0;
}

int router::assert_balance()
{
	edge_iterator it1, it2;
	double w1 = 0, w2 = 0;
	for(tie(it1, it2) = gr.in_edges(root); it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		w1 += w;
	}
	for(tie(it1, it2) = gr.out_edges(root); it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		w2 += w;
	}
	assert(fabs(w1 - w2) <= SMIN);
	return 0;
}

int router::compute_params()
{
	beta1 = beta2 = 0;
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.in_edges(root); it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		beta1 += w;
	}
	for(tie(it1, it2) = gr.out_edges(root); it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		beta1 += w;
	}
	beta1 = 1.0 / beta1 / beta1;
	beta2 = 1.0;
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

int router::build_bipartite_graph()
{
	ug.clear();
	for(int i = 0; i < u2e.size(); i++) ug.add_vertex();
	for(int i = 0; i < routes.size(); i++)
	{
		int e1 = routes[i].first;
		int e2 = routes[i].second;
		assert(e2u.find(e1) != e2u.end());
		assert(e2u.find(e2) != e2u.end());
		int s = e2u[e1];
		int t = e2u[e2];
		assert(s >= 0 && s < gr.in_degree(root));
		assert(t >= gr.in_degree(root) && t < gr.degree(root));
		ug.add_edge(s, t);
	}
	return 0;
}

int router::classify()
{
	vector<int> v = ug.assign_connected_components();

	bool b1 = true;
	bool b2 = true;
	for(int i = 1; i < gr.in_degree(root); i++)
	{
		if(v[i] != v[0]) b1 = false;
	}
	for(int i = gr.in_degree(root) + 1; i < gr.degree(root); i++)
	{
		if(v[i] != v[gr.in_degree(root)]) b2 = false;
	}
	
	if(b1 == true || b2 == true) status = INSPLITABLE;
	else status = SPLITABLE;

	return 0;
}

int router::split()
{
	if(status != SPLITABLE) return 0;

	vector< set<int> > vv = ug.compute_connected_components();

	vector<double> vw;
	for(int i = 0; i < u2e.size(); i++)
	{
		edge_descriptor e = i2e[u2e[i]];
		assert(e != null_edge);
		double w = gr.get_edge_weight(e);
		vw.push_back(w);
	}

	vector<PI> ss;
	vector<PI> tt;
	for(int i = 0; i < vv.size(); i++)
	{
		double ww = 0;
		for(set<int>::iterator it = vv[i].begin(); it != vv[i].end(); it++)
		{
			double w = vw[*it];
			if(*it >= gr.in_degree(root)) w = 0 - w;
			ww += w;
		}
		if(ww >= 0) ss.push_back(PI((int)(ww), i));
		else tt.push_back(PI((int)(0 - ww), i));
	}

	// split using subsetsum
	equation eqn1;
	equation eqn2;

	subsetsum sss(ss, tt);
	sss.solve();

	eqn1.e = sss.eqn.e;
	eqn2.e = sss.eqn.e;
	assert(eqn1.e >= 0);

	for(int i = 0; i < sss.eqn.s.size(); i++)
	{
		int k = sss.eqn.s[i];
		for(set<int>::iterator it = vv[k].begin(); it != vv[k].end(); it++)
		{
			int e = *it;
			if(e < gr.in_degree(root)) eqn1.s.push_back(u2e[e]);
			else eqn1.t.push_back(u2e[e]);
		}
	}
	for(int i = 0; i < sss.eqn.t.size(); i++)
	{
		int k = sss.eqn.t[i];
		for(set<int>::iterator it = vv[k].begin(); it != vv[k].end(); it++)
		{
			int e = *it;
			if(e < gr.in_degree(root)) eqn1.s.push_back(u2e[e]);
			else eqn1.t.push_back(u2e[e]);
		}
	}

	set<int> s1(eqn1.s.begin(), eqn1.s.end());
	set<int> s2(eqn1.t.begin(), eqn1.t.end());

	for(int i = 0; i < gr.in_degree(root); i++)
	{
		int e = u2e[i];
		if(s1.find(e) != s1.end()) continue;
		eqn2.s.push_back(e);
	}
	for(int i = gr.in_degree(root); i < gr.degree(root); i++)
	{
		int e = u2e[i];
		if(s2.find(e) != s2.end()) continue;
		eqn2.t.push_back(e);
	}

	eqns.push_back(eqn1);
	eqns.push_back(eqn2);

	delta = -2.0 * param_alpha * beta1 * eqn1.e * eqn2.e + (1.0 - param_alpha) * beta2;

	//printf("eqn.e = (%.3lf, %.3lf), alpha = %.3lf, beta = (%.3lf, %.3lf), delta = %.3lf\n", eqn1.e, eqn2.e, param_alpha, beta1, beta2, delta);

	return 0;
}

int router::decompose()
{
	if(status != INSPLITABLE) return 0;

	// build extended bipartite graph
	edge_descriptor e1 = gr.max_in_edge(root);
	edge_descriptor e2 = gr.max_out_edge(root);
	
	int k1 = -1, k2 = -1;
	for(int i = 0; i < u2e.size(); i++)
	{
		if(u2e[i] == e2i[e1]) k1 = i;
		if(u2e[i] == e2i[e2]) k2 = i;
	}
	assert(k1 != -1 && k2 != -1);

	for(int i = 0; i < gr.in_degree(root); i++)
	{
		if(ug.degree(i) >= 1) continue;
		ug.add_edge(i, k2);
	}
	for(int i = 0; i < gr.out_degree(root); i++)
	{
		int j = i + gr.in_degree(root);
		if(ug.degree(j) >= 1) continue;
		ug.add_edge(k1, j);
	}

	// run quadratic programming
	GRBEnv *env = new GRBEnv();
	GRBModel *model = new GRBModel(*env);

	// edge list of ug
	VE ve;
	edge_iterator it1, it2;
	for(tie(it1, it2) = ug.edges(); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		ve.push_back(e);
	}

	// routes~(edges in ug) weight variables
	vector<GRBVar> rvars;
	for(int i = 0; i < ve.size(); i++)
	{
		GRBVar rvar = model->addVar(1.0, GRB_INFINITY, 0, GRB_CONTINUOUS);
		rvars.push_back(rvar);
	}

	// new weights variables
	vector<GRBVar> wvars;
	for(int i = 0; i < u2e.size(); i++)
	{
		GRBVar wvar = model->addVar(1.0, GRB_INFINITY, 0, GRB_CONTINUOUS);
		wvars.push_back(wvar);
	}
	model->update();

	// expression for each edge
	vector<GRBLinExpr> exprs(u2e.size());
	for(int i = 0; i < ve.size(); i++)
	{
		edge_descriptor e = ve[i];
		int u1 = e->source();
		int u2 = e->target();
		exprs[u1] += rvars[i];
		exprs[u2] += rvars[i];
	}

	for(int i = 0; i < u2e.size(); i++)
	{
		model->addConstr(exprs[i], GRB_EQUAL, wvars[i]);
	}

	// objective 
	GRBQuadExpr obj;
	for(int i = 0; i < u2e.size(); i++)
	{
		double w = gr.get_edge_weight(i2e[u2e[i]]);
		obj += (wvars[i] - w) * (wvars[i] - w);
	}

	model->setObjective(obj, GRB_MINIMIZE);
	model->getEnv().set(GRB_IntParam_OutputFlag, 0);
	model->update();

	model->optimize();

	int f = model->get(GRB_IntAttr_Status);
	assert(f == GRB_OPTIMAL);

	if(f != GRB_OPTIMAL)
	{
		delete model;
		delete env;
		return -1;
	}

	vpi.clear();
	for(int i = 0; i < ve.size(); i++)
	{
		edge_descriptor e = ve[i];
		int s = e->source();
		int t = e->target();
		int es = u2e[s];
		int et = u2e[t];
		PI p(es, et);
		if(s > t) p = PI(et, es);
		double w = rvars[i].get(GRB_DoubleAttr_X);
		vpi.push_back(PPID(p, w));
	}

	double opt = model->get(GRB_DoubleAttr_ObjVal);
	delta = -1.0 * param_alpha * beta1 * opt + (1 - param_alpha) * beta2 * (gr.degree(root) * 1.0 - 1.0 - ve.size() * 1.0);
	//printf("opt = %.3lf, alpha = %.3lf, beta = (%.3lf, %.3lf), delta = %.3lf, degree = (%d, %lu)\n", opt, param_alpha, beta1, beta2, delta, gr.degree(root), ve.size());

	delete model;
	delete env;
	return 0;
}

int router::print() const
{
	printf("router %d, #routes = %lu, delta = %.2lf\n", root, routes.size(), delta);
	printf("in-edges = ( ");
	for(int i = 0; i < gr.in_degree(root); i++) printf("%d ", u2e[i]);
	printf("), out-edges = ( ");
	for(int i = gr.in_degree(root); i < gr.degree(root); i++) printf("%d ", u2e[i]);
	printf(")\n");

	for(int i = 0; i < routes.size(); i++)
	{
		printf("route %d (%d, %d)\n", i, routes[i].first, routes[i].second);
	}

	for(int i = 0; i < eqns.size(); i++) eqns[i].print(i);

	printf("\n");
	return 0;
}

int router::stats()
{
	vector< set<int> > vv = ug.compute_connected_components();
	
	int x1 = 0, x2 = 0;
	for(int i = 0; i < vv.size(); i++)
	{
		if(vv[i].size() <= 1) x1++;
		else x2++;
	}

	printf("vertex = %d, indegree = %d, outdegree = %d, routes = %lu, components = %lu, phased = %d, single = %d\n", 
			root, gr.in_degree(root), gr.out_degree(root), routes.size(), vv.size(), x2, x1);

	return 0;
}
