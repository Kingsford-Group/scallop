/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "router.h"
#include "config.h"
#include "util.h"
#include "subsetsum.h"

#include <iomanip>
#include <cassert>
#include <cstdio>
#include <algorithm>
#include <set>
#include <cfloat>
#include <stdint.h>

router::router(int r, splice_graph &g, MEI &ei, VE &ie)
	:root(r), gr(g), e2i(ei), i2e(ie), degree(-1), type(-1)
{
}

router::router(int r, splice_graph &g, MEI &ei, VE &ie, const MPII &mpi)
	:root(r), gr(g), e2i(ei), i2e(ie), degree(-1), type(-1)
{
	routes.clear();
	counts.clear();
	for(MPII::const_iterator it = mpi.begin(); it != mpi.end(); it++)
	{
		routes.push_back(it->first);
		counts.push_back(it->second);
	}
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

	type = rt.type;
	degree = rt.degree;
	ratio = rt.ratio;
	eqns = rt.eqns;
	pe2w = rt.pe2w;

	return (*this);
}

int router::classify()
{
	assert(gr.in_degree(root) >= 1);
	assert(gr.out_degree(root) >= 1);

	if(gr.in_degree(root) == 1 || gr.out_degree(root) == 1)
	{
		type = TRIVIAL;
		degree = gr.degree(root);
		return 0;
	}

	build_indices();
	build_bipartite_graph();
	vector< set<int> > vv = ug.compute_connected_components();

	if(routes.size() == 0)
	{
		type = SPLITTABLE_SIMPLE;
		degree = gr.degree(root) - 1;
		return 0;
	}

	if(vv.size() == 1)
	{
		type = UNSPLITTABLE_SINGLE;
		degree = ug.num_edges() - ug.num_vertices() + vv.size() + vv.size();
		return 0;
	}

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
	
	if(b1 == true || b2 == true)
	{
		type = UNSPLITTABLE_MULTIPLE;
		degree = ug.num_edges() - ug.num_vertices() + vv.size() + vv.size();
		return 0;
	}

	type = SPLITTABLE_HYPER;
	degree = vv.size() - 1;
	return 0;
}

int router::build()
{
	if(type == SPLITTABLE_SIMPLE || type == SPLITTABLE_HYPER) 
	{
		split();
	}
	if(type == UNSPLITTABLE_SINGLE || type == UNSPLITTABLE_MULTIPLE) 
	{
		thread();
	}
	return 0;
}

int router::build_indices()
{
	e2u.clear();
	u2e.clear();

	edge_iterator it1, it2;
	PEEI pei;
	for(pei = gr.in_edges(root), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int e = e2i[*it1];
		e2u.insert(PI(e, e2u.size()));
		u2e.push_back(e);
	}
	for(pei = gr.out_edges(root), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
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
	u2w.clear();
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
		edge_descriptor e = ug.add_edge(s, t);
		double w = counts[i];
		u2w.insert(PED(e, w));
	}
	return 0;
}

int router::split()
{
	eqns.clear();

	// locally smooth weights
	vector<double> vw;
	double sum1 = 0, sum2 = 0;
	for(int i = 0; i < u2e.size(); i++)
	{
		edge_descriptor e = i2e[u2e[i]];
		assert(e != null_edge);
		double w = gr.get_edge_weight(e);
		if(i < gr.in_degree(root)) sum1 += w;
		else sum2 += w;
		vw.push_back(w);
	}

	double sum = (sum1 > sum2) ? sum1 : sum2;
	double r1 = (sum1 > sum2) ? 1.0 : sum2 / sum1;
	double r2 = (sum1 < sum2) ? 1.0 : sum1 / sum2;
	
	//printf("in-degree = %d, out-degree = %d, vw.size() = %lu\n", gr.in_degree(root), gr.out_degree(root), vw.size());
	//for(int i = 0; i < vw.size(); i++) printf("debug vw[%d] = %.2lf\n", i, vw[i]);

	for(int i = 0; i < gr.in_degree(root); i++) vw[i] *= r1;
	for(int i = gr.in_degree(root); i < gr.degree(root); i++) vw[i] *= r2;

	vector< set<int> > vv = ug.compute_connected_components();

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

	// evaluate every single nontrivial component
	equation eqn0;
	eqn0.e = -1;
	for(int k = 0; k < ss.size(); k++)
	{
		set<int> &s = vv[ss[k].second];
		if(s.size() <= 1) continue;

		double r = ss[k].first * 1.0 / (sum1 * r1);
		if(eqn0.e >= 0 && r >= eqn0.e) continue;

		eqn0.clear();
		eqn0.e = r;
		assert(eqn0.e >= 0);
		for(set<int>::iterator it = s.begin(); it != s.end(); it++)
		{
			int e = *it;
			if(e < gr.in_degree(root)) eqn0.s.push_back(u2e[e]);
			else eqn0.t.push_back(u2e[e]);
		}
		assert(eqn0.s.size() >= 1);
		assert(eqn0.t.size() >= 1);
	}

	for(int k = 0; k < tt.size(); k++)
	{
		set<int> &s = vv[tt[k].second];
		if(s.size() <= 1) continue;

		double r = tt[k].first * 1.0 / (sum1 * r1);
		if(eqn0.e >= 0 && r >= eqn0.e) continue;

		eqn0.clear();
		eqn0.e = r;
		assert(eqn0.e >= 0);
		for(set<int>::iterator it = s.begin(); it != s.end(); it++)
		{
			int e = *it;
			if(e < gr.in_degree(root)) eqn0.s.push_back(u2e[e]);
			else eqn0.t.push_back(u2e[e]);
		}
		assert(eqn0.s.size() >= 1);
		assert(eqn0.t.size() >= 1);
	}

	equation eqn1;
	eqn1.e = -1;

	/*
	for(int i = 0; i < ss.size(); i++) printf("ss %d = %d:%d\n", i, ss[i].first, ss[i].second);
	for(int i = 0; i < tt.size(); i++) printf("tt %d = %d:%d\n", i, tt[i].first, tt[i].second);
	*/

	if(ss.size() >= 2 && tt.size() >= 2)
	{
		subsetsum sss(ss, tt);
		sss.solve();

		eqn1.e = sss.eqn.e;
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

		sum1 = sum2 = 0;
		for(int i = 0; i < eqn1.s.size(); i++)
		{
			int e = e2u[eqn1.s[i]];
			sum1 += vw[e];
		}
		for(int i = 0; i < eqn1.t.size(); i++)
		{
			int e = e2u[eqn1.t[i]];
			sum2 += vw[e];
		}
		eqn1.e = fabs(sum1 - sum2) / sum;
	}

	equation eqn2;
	if(eqn0.e < -0.5 && eqn1.e < -0.5) return 0;
	assert(eqn0.e >= 0 || eqn1.e >= 0);

	if(eqn1.e < -0.5) eqn2 = eqn0;
	else if(eqn0.e < -0.5) eqn2 = eqn1;
	else if(eqn0.e > eqn1.e) eqn2 = eqn1;
	else eqn2 = eqn0;

	assert(eqn2.s.size() >= 1);
	assert(eqn2.t.size() >= 1);

	set<int> s1(eqn2.s.begin(), eqn2.s.end());
	set<int> s2(eqn2.t.begin(), eqn2.t.end());

	equation eqn3;
	for(int i = 0; i < gr.in_degree(root); i++)
	{
		int e = u2e[i];
		if(s1.find(e) != s1.end()) continue;
		eqn3.s.push_back(e);
	}
	for(int i = gr.in_degree(root); i < gr.degree(root); i++)
	{
		int e = u2e[i];
		if(s2.find(e) != s2.end()) continue;
		eqn3.t.push_back(e);
	}

	if(eqn3.s.size() <= 0 || eqn3.t.size() <= 0) return 0;

	ratio = eqn2.e;
	eqn3.e = ratio;

	eqns.push_back(eqn2);
	eqns.push_back(eqn3);

	return 0;
}

int router::thread()
{
	pe2w.clear();
	vector<int> v1;
	for(int k = 0; k < u2e.size(); k++)
	{
		int e = u2e[k];
		if(ug.degree(k) == 0) v1.push_back(k);
	}

	vector<double> vw = compute_balanced_weights();
	double weight_sum = 0;
	for(int k = 0; k < vw.size(); k++) weight_sum += vw[k];

	bool b;
	while(true)
	{
		b = thread_leaf(vw);
		if(b == true) continue;

		b = thread_turn(vw);
		if(b == false) break;
	}

	assert(ug.num_edges() == 0);

	int n = gr.in_degree(root);
	for(int i = 0; i < v1.size(); i++)
	{
		int k = v1[i];
		if(k < n) thread_isolate1(k, vw);
		else thread_isolate2(k, vw);
	}
	
	double weight_remain = 0;
	for(int k = 0; k < vw.size(); k++)
	{
		//printf("weight remain for edge %d = %.2lf, sum = %.2lf\n", u2e[k], vw[k], weight_sum);
		if(vw[k] <= 0) continue;
		weight_remain += vw[k];
	}

	ratio = weight_remain / weight_sum;

	for(MPID::iterator it = pe2w.begin(); it != pe2w.end(); it++)
	{
		if(it->second < 1.0) it->second = 1.0;
	}
	return 0;
}

bool router::thread_leaf(vector<double> &vw)
{
	PEEI pei = ug.edges();
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		int s = (*it)->source();
		int t = (*it)->target();

		if(s >= t)
		{
			int a = s;
			s = t;
			t = a;
		}

		if(vw[s] < -0.5) continue;
		if(vw[t] < -0.5) continue;

		if(ug.degree(s) == 1 && vw[s] <= vw[t])
		{
			PPID pw(PI(u2e[s], u2e[t]), vw[s]);
			pe2w.insert(pw);
			ug.clear_vertex(s);
			vw[t] -= vw[s];
			vw[s] = -1;
			return true;
		}
		if(ug.degree(t) == 1 && vw[t] <= vw[s])
		{
			PPID pw(PI(u2e[s], u2e[t]), vw[t]);
			pe2w.insert(pw);
			ug.clear_vertex(t);
			vw[s] -= vw[t];
			vw[t] = -1;
			return true;
		}
	}
	return false;
}

bool router::thread_turn(vector<double> &vw)
{
	int x = -1;
	for(int k = 0; k < vw.size(); k++)
	{
		if(vw[k] < -0.5) continue;
		if(ug.degree(k) <= 1) continue;
		if(x != -1 && vw[k] > vw[x]) continue;
		x = k;
	}

	if(x == -1) return false;

	double sum = 0;
	PEEI pei = ug.out_edges(x);
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		int s = (*it)->source();
		int t = (*it)->target();
		assert(s == x);
		sum += u2w[*it];
		assert(vw[t] >= vw[x]);
	}

	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		int t = (*it)->target();
		double w = vw[x] * u2w[*it] / sum;
		//;if(u2w[*it] == 1) w = 1;	// set to 1 for those with only 1 read supported
		PI p = (x < t) ? PI(u2e[x], u2e[t]) : PI(u2e[t], u2e[x]);
		PPID pw(p, w);
		pe2w.insert(pw);
		vw[t] -= w;
	}

	vw[x] = -1;
	ug.clear_vertex(x);
	return true;
}

int router::thread_isolate1(int k, vector<double> &vw)
{
	int x = gr.in_degree(root);
	for(int i = gr.in_degree(root) + 1; i < u2e.size(); i++)
	{
		if(vw[i] < vw[x]) continue;
		x = i;
	}
	assert(x != -1);
	double w = vw[x] < vw[k] ? vw[x] : vw[k];
	vw[x] -= w;
	vw[k] -= w;
	PPID pw(PI(u2e[k], u2e[x]), w);
	pe2w.insert(pw);
	return 0;
}

int router::thread_isolate2(int k, vector<double> &vw)
{
	int x = 0;
	for(int i = 0; i < gr.in_degree(root); i++)
	{
		if(vw[i] < vw[x]) continue;
		x = i;
	}
	assert(x != -1);
	double w = vw[x] < vw[k] ? vw[x] : vw[k];
	vw[x] -= w;
	vw[k] -= w;
	PPID pw(PI(u2e[x], u2e[k]), w);
	pe2w.insert(pw);
	return 0;
}

vector<double> router::compute_balanced_weights()
{
	vector<double> vw;
	double sum1 = 0, sum2 = 0;
	for(int i = 0; i < u2e.size(); i++)
	{
		edge_descriptor e = i2e[u2e[i]];
		assert(e != null_edge);
		double w = gr.get_edge_weight(e);
		if(i < gr.in_degree(root)) sum1 += w;
		else sum2 += w;
		vw.push_back(w);
	}
	double r1 = sqrt(sum2 / sum1);
	double r2 = sqrt(sum1 / sum2);
	for(int i = 0; i < gr.in_degree(root); i++) vw[i] *= r1;
	for(int i = gr.in_degree(root); i < gr.degree(root); i++) vw[i] *= r2;
	
	return vw;
}

int router::print()
{
	vector<double> vw = compute_balanced_weights();

	printf("router %d, #routes = %lu, type = %d, degree = %d, ratio = %.2lf\n", root, routes.size(), type, degree, ratio);
	printf("in-edges = ( ");
	for(int i = 0; i < gr.in_degree(root); i++) printf("%d ", u2e[i]);
	printf("), weights = ( ");
	for(int i = 0; i < gr.in_degree(root); i++) printf("%.1lf,%.1lf ", gr.get_edge_weight(i2e[u2e[i]]), vw[i]);
	printf(")\n");

	printf("out-edges = ( ");
	for(int i = gr.in_degree(root); i < gr.degree(root); i++) printf("%d ", u2e[i]);
	printf("), weights = ( ");
	for(int i = gr.in_degree(root); i < gr.degree(root); i++) printf("%.1lf,%.1lf ", gr.get_edge_weight(i2e[u2e[i]]), vw[i]);
	printf(")\n");

	for(int i = 0; i < routes.size(); i++)
	{
		printf("route %d (%d, %d), count = %d\n", i, routes[i].first, routes[i].second, counts[i]);
	}

	for(int i = 0; i < eqns.size(); i++) eqns[i].print(i);

	for(MPID::const_iterator it = pe2w.begin(); it != pe2w.end(); it++)
	{
		printf("decompose: (%d, %d), w = %.2lf\n", it->first.first, it->first.second, it->second);
	}

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

bool compare_edge_weight(const PED &x, const PED &y)
{
	if(x.second > y.second) return true;
	else return false;
}
