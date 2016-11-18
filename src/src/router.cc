#include "router.h"
#include "util.h"
#include "gurobi_c++.h"
#include "smoother.h"
#include "subsetsum.h"
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

	bw = rt.bw;
	bratio = rt.bratio;

	ratio = rt.ratio;
	eqns = rt.eqns;

	status = rt.status;

	return (*this);
}

int router::build()
{
	assert(gr.in_degree(root) >= 1);
	assert(gr.out_degree(root) >= 1);

	eqns.clear();
	ratio = -1;
	status = -1;

	build_indices();

	if(gr.in_degree(root) == 1 || gr.out_degree(root) == 1)
	{
		add_single_equation();
		status = 0;
		return 0;
	}

	build_bipartite_graph();
	vector< set<int> > vv = ug.compute_connected_components();

	if(vv.size() == 1)
	{
		if(ug.num_edges() == ug.num_vertices() - 1) status = 1;
		else if(ug.num_edges() >= u2e.size()) status = 2;
		else assert(false);
		return 0;
	}

	bool b1 = true;
	bool b2 = true;
	vector<int> v = ug.assign_connected_components();
	for(int i = 1; i < gr.in_degree(root); i++)
	{
		if(v[i] != v[0]) b1 = false;
	}
	for(int i = gr.in_degree(root) + 1; i < gr.degree(root); i++)
	{
		if(v[i] != v[gr.in_degree(root)]) b2 = false;
	}
	
	if(b1 == true || b2 == true)
	{
		status = 5;
		return 0;
	}

	if(vv.size() == 2) status = 3;
	else status = 4;

	split();

	assert(eqns.size() == 2);

	return 0;
}

int router::add_single_equation()
{
	equation eqn;
	double sum1 = 0, sum2 = 0;
	for(int i = 0; i < gr.in_degree(root); i++)
	{
		eqn.s.push_back(u2e[i]);
		sum1 += gr.get_edge_weight(i2e[u2e[i]]);
	}
	for(int i = gr.in_degree(root); i < gr.degree(root); i++)
	{
		eqn.t.push_back(u2e[i]);
		sum2 += gr.get_edge_weight(i2e[u2e[i]]);
	}
	assert(eqn.s.size() >= 1);
	assert(eqn.t.size() >= 1);
	
	ratio = fabs(sum1 - sum2) / (sum1 + sum2);
	eqn.e = ratio;

	eqns.push_back(eqn);
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

int router::build_balanced_weights()
{
	set<int> fb;
	return build_balanced_weights(fb);
}

int router::build_balanced_weights(const set<int> &fb)
{
	bw.assign(u2e.size(), 0);
	double sum1 = 0, sum2 = 0;
	for(int i = 0; i < u2e.size(); i++)
	{
		if(fb.find(i) != fb.end()) continue;
		edge_descriptor e = i2e[u2e[i]];
		assert(e != null_edge);
		double w = gr.get_edge_weight(e);
		bw[i] = w;
		if(i < gr.in_degree(root)) sum1 += w;
		else sum2 += w;
	}

	if(sum1 >= sum2) bratio = sqrt(sum1 / sum2);
	else bratio = sqrt(sum2 / sum1);

	for(int i = 0; i < gr.in_degree(root); i++)
	{
		if(sum1 >= sum2) bw[i] /= bratio;
		else bw[i] *= bratio;
	}

	for(int i = gr.in_degree(root); i < gr.degree(root); i++)
	{
		if(sum2 >= sum1) bw[i] /= bratio;
		else bw[i] *= bratio;
	}
	return 0;
}

int router::split()
{
	vector< set<int> > vv = ug.compute_connected_components();

	build_balanced_weights();

	vector<PI> ss;
	vector<PI> tt;
	for(int i = 0; i < vv.size(); i++)
	{
		double ww = 0;
		for(set<int>::iterator it = vv[i].begin(); it != vv[i].end(); it++)
		{
			double w = bw[*it];
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

	ratio = (eqn1.e * 0.5 + eqn2.e * 0.5);
	//ratio = (eqn1.e * 0.5 + eqn2.e * 0.5) * bratio;

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

int router::print() const
{
	printf("router %d, #routes = %lu, ratio = %.2lf\n", root, routes.size(), ratio);
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
