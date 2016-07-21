#include "scallop2.h"
#include "subsetsum.h"
#include "config.h"
#include "smoother.h"
#include "nested_graph.h"

#include <cstdio>
#include <iostream>
#include <cfloat>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/adjacency_list.hpp>

using namespace boost;

scallop2::scallop2()
{}

scallop2::scallop2(const string &s, splice_graph &g)
	: name(s), gr(g)
{
	round = 0;
	//assert(gr.check_fully_connected() == true);
}

scallop2::~scallop2()
{}

int scallop2::clear()
{
	name = "";
	gr.clear();
	e2i.clear();
	i2e.clear();
	mev.clear();
	round = -1;
	paths.clear();
	return 0;
}

int scallop2::save(scallop2 &sc)
{
	sc.clear();

	sc.name = name;
	MEE x2y, y2x;
	sc.gr.copy(gr, x2y, y2x);

	for(MEI::iterator it = e2i.begin(); it != e2i.end(); it++)
	{
		edge_descriptor e = it->first;
		int k = it->second;
		assert(k >= 0 && k < i2e.size());
		if(i2e[k] == null_edge)
		{
			assert(x2y.find(e) == x2y.end());
			sc.e2i.insert(PEI(e, k));
		}
		else
		{
			assert(x2y.find(e) != x2y.end());
			sc.e2i.insert(PEI(x2y[e], k));
		}
	}

	for(int i = 0; i < i2e.size(); i++)
	{
		if(i2e[i] == null_edge)
		{
			sc.i2e.push_back(null_edge);
		}
		else
		{
			assert(x2y.find(i2e[i]) != x2y.end());
			sc.i2e.push_back(x2y[i2e[i]]);
		}
	}

	for(MEV::iterator it = mev.begin(); it != mev.end(); it++)
	{
		edge_descriptor e = it->first;
		vector<int> v = it->second;
		if(x2y.find(e) == x2y.end())
		{
			if(e2i.find(e) == e2i.end()) continue;
			int k = e2i[e];
			assert(i2e[k] == null_edge);
			assert(sc.i2e[k] == null_edge);
			sc.mev.insert(PEV(e, v));
		}
		else
		{
			int k = e2i[e];
			assert(i2e[k] == e);
			assert(sc.e2i[x2y[e]] == k);
			assert(sc.i2e[k] == x2y[e]);
			sc.mev.insert(PEV(x2y[e], v));
		}
	}

	sc.round = round;
	sc.paths = paths;

	return 0;
}

int scallop2::load(scallop2 &sc)
{
	clear();

	name = sc.name;

	MEE x2y, y2x;
	gr.copy(sc.gr, x2y, y2x);

	for(MEI::iterator it = sc.e2i.begin(); it != sc.e2i.end(); it++)
	{
		edge_descriptor e = it->first;
		int k = it->second;
		assert(k >= 0 && k < sc.i2e.size());
		if(sc.i2e[k] == null_edge)
		{
			assert(x2y.find(e) == x2y.end());
			e2i.insert(PEI(e, k));
		}
		else
		{
			assert(x2y.find(e) != x2y.end());
			e2i.insert(PEI(x2y[e], k));
		}
	}

	for(int i = 0; i < sc.i2e.size(); i++)
	{
		if(sc.i2e[i] == null_edge)
		{
			i2e.push_back(null_edge);
		}
		else
		{
			assert(x2y.find(sc.i2e[i]) != x2y.end());
			i2e.push_back(x2y[sc.i2e[i]]);
		}
	}

	for(MEV::iterator it = sc.mev.begin(); it != sc.mev.end(); it++)
	{
		edge_descriptor e = it->first;
		vector<int> v = it->second;
		if(x2y.find(e) == x2y.end())
		{
			if(sc.e2i.find(e) == sc.e2i.end()) continue;
			int k = sc.e2i[e];
			assert(sc.i2e[k] == null_edge);
			assert(i2e[k] == null_edge);
			mev.insert(PEV(e, v));
		}
		else
		{
			int k = sc.e2i[e];
			assert(sc.i2e[k] == e);
			assert(e2i[x2y[e]] == k);
			assert(i2e[k] == x2y[e]);
			mev.insert(PEV(x2y[e], v));
		}
	}

	round = sc.round;
	paths = sc.paths;

	return 0;
}

int scallop2::assemble()
{
	int c = classify();
	//if(c == TRIVIAL) return 0;

	if(algo == "core") return assemble1();
	if(algo == "full") return assemble2();
	if(algo == "greedy") return greedy();
	return 0;
}

int scallop2::classify()
{
	assert(gr.num_vertices() >= 2);
	if(gr.num_vertices() == 2) return TRIVIAL;

	string s;	

	long p0 = gr.compute_num_paths();
	long p1 = gr.num_edges() - gr.num_vertices() + 2;
	for(int i = 0; i < gr.num_vertices(); i++) 
	{
		if(gr.degree(i) == 0) p1++;
	}

	printf("vertices = %lu, edges = %lu, p0 = %ld, p1 = %ld\n", gr.num_vertices(), gr.num_edges(), p0, p1);

	assert(p0 >= p1);

	bool b = (p0 == p1) ? true : false;

	printf("\nprocess %s %s\n", name.c_str(), b ? "TRIVIAL" : "NORMAL");

	if(p0 == p1) return TRIVIAL;
	else return NORMAL;
}

int scallop2::assemble0()
{
	if(output_tex_files == true) gr.draw(name + "." + tostring(round++) + ".tex");

	gr.get_edge_indices(i2e, e2i);
	init_super_edges();

	print();

	while(true)
	{
		bool b = false;

		b = join_trivial_edges();
		if(b == true) print();
		if(b == true) continue;

		b = join_trivial_vertices();
		if(b == true) print();
		if(b == true) continue;
		
		break;
	}

	smooth();

	return 0;
}

int scallop2::assemble1()
{
	assemble0();
	int f = iterate();

	collect_existing_st_paths();

	printf("%s core solution %lu paths, iteration = %d\n", name.c_str(), paths.size(), f);

	return 0;
}

int scallop2::assemble2()
{
	assemble0();
	int f = iterate();

	greedy_decompose();
	assert(gr.num_edges() == 0);

	printf("%s full solution %lu paths, iteration = %d\n", name.c_str(), paths.size(), f);

	return 0;
}

int scallop2::greedy()
{
	assemble0();

	greedy_decompose();
	assert(gr.num_edges() == 0);

	printf("%s greedy solution %lu paths\n", name.c_str(), paths.size());

	return 0;
}

int scallop2::iterate()
{
	while(true)
	{
		bool b = false;

		b = decompose_trivial_internal_vertices();
		if(b == true) print();
		if(b == true) continue;

		b = decompose_with_equations(1);
		if(b == true) print();
		if(b == true) continue;

		break;

		b = decompose_with_equations(2);
		if(b == true) print();
		if(b == true) continue;

		break;
	}

	bool flag = true;
	while(true)
	{
		bool b = false;
		b = decompose_trivial_external_vertices();
		if(b == true) flag = true;
		if(b == true) continue;

		b = decompose_trivial_internal_vertices();
		if(b == true) flag = true;
		if(b == true) continue;

		break;
	}

	if(flag == true) print();

	return 0;
}

bool scallop2::decompose_with_equations(int level)
{
	vector<equation> eqns;

	//gr.round_weights();

	bool b = false;
	if(level == 1) b = identify_equations1(eqns);
	if(level == 2) b = identify_equations2(eqns);

	if(eqns.size() == 0) return false;

	sort(eqns.begin(), eqns.end(), equation_cmp1);

	printf("candidate equations: %lu\n", eqns.size());
	/*
	for(int i = 0; i < eqns.size(); i++)
	{
		eqns[i].print(i);
	}
	*/

	equation eqn(0);
	for(int i = 0; i < eqns.size(); i++)
	{
		eqn = eqns[i];

		if(verify_equation_mergable(eqn) == false) continue;
		if(verify_equation_nontrivial(eqn) == false) continue;

		eqns[i].print(i);
		scallop2 sc;
		save(sc);
		bool b = smooth(eqn);
		if(b == true) resolve_equation(eqn);
		load(sc);

		if(eqn.f == 2 && eqn.d == 0) break;
	}

	if(eqn.f != 2 || eqn.d != 0) return false;

	printf("smooth with equation\n");
	eqn.print(99);
	smooth(eqn);
	print();

	printf("resolve with equation\n");
	eqn.print(99);
	resolve_equation(eqn);

	return true;
}

int scallop2::init_super_edges()
{
	mev.clear();
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		vector<int> v;
		int s = (*it1)->source();
		v.push_back(s);
		mev.insert(PEV(*it1, v));
	}
	return 0;
}

int scallop2::split_merge_path(const VE &p, double wx, vector<int> &vv)
{
	vector<int> v;
	for(int i = 0; i < p.size(); i++)
	{
		assert(p[i] != null_edge);
		assert(e2i.find(p[i]) != e2i.end());
		v.push_back(e2i[p[i]]);
	}
	return split_merge_path(v, wx, vv);
}

int scallop2::split_merge_path(const vector<int> &p, double ww, vector<int> &vv)
{
	vv.clear();
	if(p.size() == 0) return -1;
	int ee = p[0];
	int x = split_edge(p[0], ww);
	vv.push_back(x);
	for(int i = 1; i < p.size(); i++)
	{
		x = split_edge(p[i], ww);
		vv.push_back(x);
		ee = merge_adjacent_equal_edges(ee, p[i]);
	}
	return ee;
}

int scallop2::merge_adjacent_equal_edges(int x, int y)
{
	if(i2e[x] == null_edge) return -1;
	if(i2e[y] == null_edge) return -1;

	edge_descriptor xx = i2e[x];
	edge_descriptor yy = i2e[y];

	int xs = (xx)->source();
	int xt = (xx)->target();
	int ys = (yy)->source();
	int yt = (yy)->target();

	if(xt != ys && yt != xs) return -1;
	if(yt == xs) return merge_adjacent_equal_edges(y, x);
	
	assert(xt == ys);

	edge_descriptor p = gr.add_edge(xs, yt);

	int n = i2e.size();
	i2e.push_back(p);
	assert(e2i.find(p) == e2i.end());
	e2i.insert(PEI(p, n));

	double wx0 = gr.get_edge_weight(xx);
	double wy0 = gr.get_edge_weight(yy);
	double wx1 = gr.get_edge_stddev(xx);
	double wy1 = gr.get_edge_stddev(yy);

	assert(fabs(wx0 - wy0) <= SMIN);

	gr.set_edge_weight(p, wx0);
	gr.set_edge_stddev(p, wx1 + wy1);

	vector<int> v = mev[xx];
	v.insert(v.end(), mev[yy].begin(), mev[yy].end());

	if(mev.find(p) != mev.end()) mev[p] = v;
	else mev.insert(PEV(p, v));

	assert(i2e[n] == p);
	assert(e2i.find(p) != e2i.end());
	assert(e2i[p] == n);
	assert(e2i[i2e[n]] == n);

	e2i.erase(xx);
	e2i.erase(yy);
	i2e[x] = null_edge;
	i2e[y] = null_edge;
	gr.remove_edge(xx);
	gr.remove_edge(yy);

	return n;
}

int scallop2::merge_adjacent_edges(int x, int y)
{
	if(i2e[x] == null_edge) return -1;
	if(i2e[y] == null_edge) return -1;

	edge_descriptor xx = i2e[x];
	edge_descriptor yy = i2e[y];

	double wx = gr.get_edge_weight(xx);
	double wy = gr.get_edge_weight(yy);
	double ww = (wx <= wy) ? wx : wy;

	split_edge(x, ww);
	split_edge(y, ww);

	return merge_adjacent_equal_edges(x, y);
}

int scallop2::split_edge(int ei, double w)
{
	assert(i2e[ei] != null_edge);
	edge_descriptor ee = i2e[ei];

	double ww = gr.get_edge_weight(ee);
	double dd = gr.get_edge_stddev(ee);

	if(fabs(ww - w) <= SMIN) return ei;
	assert(ww >= w + SMIN);

	int s = ee->source();
	int t = ee->target();

	edge_descriptor p2 = gr.add_edge(s, t);

	gr.set_edge_weight(ee, w);
	gr.set_edge_stddev(ee, dd);
	gr.set_edge_weight(p2, ww - w);
	gr.set_edge_stddev(p2, dd);

	if(mev.find(p2) != mev.end()) mev[p2] = mev[ee];
	else mev.insert(PEV(p2, mev[ee]));

	int n = i2e.size();
	i2e.push_back(p2);
	e2i.insert(PEI(p2, n));

	return n;
}

bool scallop2::verify_equation_mergable(equation &eqn)
{
	if(eqn.s.size() == 0) return false; 
	if(eqn.t.size() == 0) return false;

	SE xx, yy;
	for(int i = 0; i < eqn.s.size(); i++)
	{
		edge_descriptor x = i2e[eqn.s[i]];
		SE se;
		gr.bfs(x->target(), se);
		xx.insert(se.begin(), se.end());
		se.clear();
		gr.bfs_reverse(x->source(), se);
		xx.insert(se.begin(), se.end());
	}
	for(int i = 0; i < eqn.t.size(); i++)
	{
		edge_descriptor x = i2e[eqn.t[i]];
		SE se;
		gr.bfs(x->target(), se);
		yy.insert(se.begin(), se.end());
		se.clear();
		gr.bfs_reverse(x->source(), se);
		yy.insert(se.begin(), se.end());
	}

	for(int i = 0; i < eqn.s.size(); i++)
	{
		edge_descriptor x = i2e[eqn.s[i]];
		if(yy.find(x) == yy.end()) return false;
	}

	for(int i = 0; i < eqn.t.size(); i++)
	{
		edge_descriptor x = i2e[eqn.t[i]];
		if(xx.find(x) == xx.end()) return false;
	}
	return true;
}

bool scallop2::verify_equation_nontrivial(equation &eqn)
{
	set<edge_descriptor> fbs;
	set<edge_descriptor> fbt;
	vector<int> vsx;
	vector<int> vsy;
	vector<int> vtx;
	vector<int> vty;
	for(int i = 0; i < eqn.s.size(); i++)
	{
		edge_descriptor &e = i2e[eqn.s[i]];
		assert(e != null_edge);
		fbs.insert(e);
		vsx.push_back(e->source());
		vsy.push_back(e->target());
	}
	for(int i = 0; i < eqn.t.size(); i++)
	{
		edge_descriptor &e = i2e[eqn.t[i]];
		assert(e != null_edge);
		fbt.insert(e);
		vtx.push_back(e->source());
		vty.push_back(e->target());
	}

	int n = gr.num_vertices() - 1;
	bool b1 = gr.bfs(vsy, n, fbt);
	bool b2 = gr.bfs_reverse(vtx, 0, fbs);
	if(b1 == false || b2 == false) return false;

	b1 = gr.bfs(vty, n, fbs);
	b2 = gr.bfs_reverse(vsx, 0, fbt);
	if(b1 == false || b2 == false) return false;

	return true;
}

int scallop2::identify_equation(const vector<int> &subs, vector<equation> &eqns)
{
	if(subs.size() == 0) return 0;

	int e = subs[0];
	int s = i2e[e]->source();
	int t = i2e[e]->target();
	double w = gr.get_edge_weight(i2e[e]);

	SE ff, bb;
	gr.bfs(t, ff);
	gr.bfs_reverse(s, bb);
	ff.insert(bb.begin(), bb.end());

	VE vv(gr.num_edges());

	for(int i = 1; i < subs.size(); i++)
	{
		e = subs[i];
		s = i2e[e]->source();
		t = i2e[e]->target();
		w += gr.get_edge_weight(i2e[e]);

		SE f, b;
		gr.bfs(t, f);
		gr.bfs_reverse(s, b);
		f.insert(b.begin(), b.end());

		VE::iterator it = set_union(f.begin(), f.end(), ff.begin(), ff.end(), vv.begin());
		ff = SE(vv.begin(), it);
	}

	int sw = (int)(w);

	set<int> ss(subs.begin(), subs.end());
	vector<PI> xi;
	for(SE::iterator it = ff.begin(); it != ff.end(); it++)
	{
		edge_descriptor e = *it;
		if(ss.find(e2i[e]) != ss.end()) continue;
		int s = e->source();
		int t = e->target();
		if(s == 0 || t == gr.num_vertices() - 1) continue;
		int ww = (int)(gr.get_edge_weight(*it));
		if(ww <= 0) continue;
		if(ww > sw) continue;
		//if(ww * 1.0 > sw * 1.0 * (1.0 + max_equation_error_ratio)) continue;
		xi.push_back(PI(ww, e2i[*it]));
	}

	if(xi.size() == 0) return 0;

	sort(xi.begin(), xi.end());
	xi.push_back(PI(sw, e));

	vector<int> v;
	for(int i = 0; i < xi.size(); i++)
	{
		v.push_back(xi[i].first);
	}

	subsetsum sss(v);
	sss.solve();

	if(sss.subsets.size() == 0) return 0;

	for(int i = 0; i < sss.subsets.size(); i++)
	{
		vector<int> subset = sss.subsets[i];
		int opt = sss.opts[i];

		vector<int> subt;
		double wt = 0;
		for(int j = 0; j < subset.size(); j++)
		{
			int k = subset[j];
			int e = xi[k].second;
			subt.push_back(e);
			wt += gr.get_edge_weight(i2e[e]);
		}

		double err = fabs(wt - w);

		if(err * 1.0 > 1.0 * sw * max_equation_error_ratio) continue;

		equation eqn(err);
		eqn.s = subs;
		eqn.t = subt;
		eqns.push_back(eqn);

		//return 0;		// only consider one equation
	}
	return 0;
}

int scallop2::identify_equations1(vector<equation> &eqns)
{
	for(int i = 0; i < i2e.size(); i++)
	{
		edge_descriptor &e = i2e[i];
		if(e == null_edge) continue;
		if(e->source() == 0) continue;
		if(e->target() == gr.num_vertices() - 1) continue;

		vector<int> s;
		s.push_back(i);

		identify_equation(s, eqns);
	}
	return 0;

}

int scallop2::identify_equations2(vector<equation> &eqns)
{
	vector<int> p;
	for(int i = 0; i < i2e.size(); i++)
	{
		edge_descriptor &e = i2e[i];
		if(e == null_edge) continue;
		if(e->source() == 0) continue;
		if(e->target() == gr.num_vertices() - 1) continue;
		int w = (int)(gr.get_edge_weight(e));
		if(w <= 0) continue;
		p.push_back(i);
	}

	if(p.size() == 0) return 0;

	for(int i = 0; i < p.size(); i++)
	{
		for(int j = i + 1; j < p.size(); j++)
		{
			vector<int> s;
			s.push_back(p[i]);
			s.push_back(p[j]);
			
			identify_equation(s, eqns);
		}
	}

	return 0;
}

bool scallop2::resolve_vertex_with_equations(vector<equation> &eqns)
{
	for(int i = 0; i < eqns.size(); i++)
	{
		equation &eqn = eqns[i];
		resolve_vertex_with_equation(eqn);
		if(eqn.f == 3) return true;
	}
	return false;
}

bool scallop2::resolve_vertex_with_equation(equation &eqn)
{
	resolve_vertex_with_equation1(eqn);
	resolve_vertex_with_equation2(eqn);
	if(eqn.f == 3) return true;

	vector<int> t = eqn.s;
	eqn.s = eqn.t;
	eqn.t = t;
	resolve_vertex_with_equation1(eqn);
	resolve_vertex_with_equation2(eqn);
	if(eqn.f == 3) return true;
	return false;
}

bool scallop2::resolve_vertex_with_equation1(equation &eqn)
{
	if(eqn.s.size() == 0) return false;	
	
	SE sex;
	edge_descriptor ex = i2e[eqn.s[0]];
	sex.insert(ex);
	int xt = ex->target();
	for(int i = 1; i < eqn.s.size(); i++)
	{
		edge_descriptor ex = i2e[eqn.s[i]];
		sex.insert(ex);
		if(ex->target() != xt) return false;
	}

	for(int i = 0; i < eqn.t.size(); i++)
	{
		edge_descriptor ey = i2e[eqn.t[i]];
		int ys = ey->source();
		if(gr.check_path(xt, ys) == false) return false;
	}

	assert(gr.in_degree(xt) >= 1);
	if(gr.in_degree(xt) != 1 + eqn.s.size()) return false;

	edge_descriptor ex2 = null_edge;
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.in_edges(xt); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		if(sex.find(e) != sex.end()) continue;
		ex2 = e;
	}
	assert(ex2 != null_edge);
	double w2 = gr.get_edge_weight(ex2);

	SE sey;
	for(int i = 0; i < eqn.t.size(); i++)
	{
		edge_descriptor ey = i2e[eqn.t[i]];
		sey.insert(ey);
		int ys = ey->source();
		SE se;
		gr.bfs_reverse(ys, se);
		sey.insert(se.begin(), se.end());
	}

	VE ve;
	double w1 = 0;
	for(tie(it1, it2) = gr.out_edges(xt); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		if(sey.find(e) != sey.end()) continue;
		ve.push_back(e);
		w1 += gr.get_edge_weight(e);
	}

	if(ve.size() <= 0) return false;
	if(ve.size() >= 1 && w1 > w2 + SMIN) return false;

	int ei = e2i[ex2];
	for(int i = 0; i < ve.size(); i++)
	{
		edge_descriptor e = ve[i];
		double ww = gr.get_edge_weight(e);
		int ee = split_edge(ei, ww);

		printf("resolve vertex with equation: (%d, %d)\n", ei, e2i[e]);

		merge_adjacent_equal_edges(ei, e2i[e]);
		ei = ee;
	}

	eqn.f = 3;
	return true;
}

bool scallop2::resolve_vertex_with_equation2(equation &eqn)
{
	if(eqn.t.size() == 0) return false;	

	SE sey;
	edge_descriptor ey = i2e[eqn.t[0]];
	sey.insert(ey);
	int ys = ey->source();
	for(int i = 1; i < eqn.t.size(); i++)
	{
		edge_descriptor ey = i2e[eqn.t[i]];
		sey.insert(ey);
		if(ey->source() != ys) return false;
	}

	for(int i = 0; i < eqn.s.size(); i++)
	{
		edge_descriptor ex = i2e[eqn.s[i]];
		int xt = ex->target();
		if(gr.check_path(xt, ys) == false) return false;
	}

	assert(gr.out_degree(ys) >= 1);
	if(gr.out_degree(ys) != 1 + eqn.t.size()) return false;

	edge_iterator it1, it2;
	edge_descriptor ey2 = null_edge;
	for(tie(it1, it2) = gr.out_edges(ys); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		if(sey.find(e) != sey.end()) continue;
		ey2 = e;
	}
	assert(ey2 != null_edge);
	double w2 = gr.get_edge_weight(ey2);

	SE sex;
	for(int i = 0; i < eqn.s.size(); i++)
	{
		edge_descriptor ex = i2e[eqn.s[i]];
		sex.insert(ex);
		int xt = ex->target();
		SE se;	
		gr.bfs(xt, se);
		sex.insert(se.begin(), se.end());
	}
	VE ve;
	double w1 = 0;
	for(tie(it1, it2) = gr.in_edges(ys); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		if(sex.find(e) != sex.end()) continue;
		ve.push_back(e);
		w1 += gr.get_edge_weight(e);
	}

	if(ve.size() <= 0) return false;
	if(ve.size() >= 1 && w1 > w2 + SMIN) return false;

	int ei = e2i[ey2];
	for(int i = 0; i < ve.size(); i++)
	{
		edge_descriptor e = ve[i];
		double ww = gr.get_edge_weight(e);
		int ee = split_edge(ei, ww);

		printf("resolve vertex with equation: (%d, %d)\n", ei, e2i[e]);

		merge_adjacent_equal_edges(ei, e2i[e]);
		ei = ee;
	}

	eqn.f = 3;
	return true;
}

int scallop2::smooth()
{
	smoother sm(gr);
	sm.solve();
	return 0;
}

bool scallop2::smooth(equation &eqn)
{
	VE vx, vy;
	for(int i = 0; i < eqn.s.size(); i++) vx.push_back(i2e[eqn.s[i]]);
	for(int i = 0; i < eqn.t.size(); i++) vy.push_back(i2e[eqn.t[i]]);

	smoother sm(gr);
	sm.add_equation(vx, vy);
	bool f = sm.solve();
	if(f == 0) return true;
	else return false;
}

int scallop2::resolve_equation(equation &eqn)
{
	vector<int> s = eqn.s;
	vector<int> t = eqn.t;
	eqn.a = eqn.d = 0;
	eqn.f = resolve_equation(s, t, eqn.a, eqn.d);
	return 0;
}

int scallop2::resolve_equation(vector<int> &s, vector<int> &t, int &ma, int &md)
{
	if(s.size() == 0 && t.size() == 0) return 2;

	assert(s.size() >= 1);
	assert(t.size() >= 1);

	for(int i = 0; i < s.size(); i++)
	{
		for(int j = 0; j < t.size(); j++)
		{
			int x = s[i];
			int y = t[j];
			assert(i2e[x] != null_edge);
			assert(i2e[y] != null_edge);

			vector<PI> p;
			if(check_adjacent_mergable(x, y, p) == false) continue;

			build_adjacent_edges(p);

			double wx = gr.get_edge_weight(i2e[x]);
			double wy = gr.get_edge_weight(i2e[y]);
			double ww = (wx <= wy) ? wx : wy;

			vector<int> v;
			v.push_back(x);
			v.push_back(y);

			vector<int> vv;
			split_merge_path(v, ww, vv);

			ma++;

			//// printf("merge (adjacent) edge pair (%d, %d)\n", x, y);

			if(i2e[vv[0]] == null_edge) assert(vv[0] == x);
			if(i2e[vv[1]] == null_edge) assert(vv[1] == y);

			if(vv[0] == x) s.erase(s.begin() + i);
			else s[i] = vv[0];

			if(vv[1] == y) t.erase(t.begin() + j);
			else t[j] = vv[1];

			int f = resolve_equation(s, t, ma, md);
			if(f == 2) return 2;
			else return 1;
		}
	}

	set<int> ss;
	for(int i = 0; i < s.size(); i++)
	{
		assert(ss.find(s[i]) == ss.end());
		ss.insert(s[i]);
	}
	for(int i = 0; i < t.size(); i++)
	{
		assert(ss.find(t[i]) == ss.end());
		ss.insert(t[i]);
	}

	for(int i = 0; i < s.size(); i++)
	{
		for(int j = 0; j < t.size(); j++)
		{
			int x = s[i];
			int y = t[j];
			assert(i2e[x] != null_edge);
			assert(i2e[y] != null_edge);

			double wx = gr.get_edge_weight(i2e[x]);
			double wy = gr.get_edge_weight(i2e[y]);
			double ww = (wx <= wy) ? wx : wy;

			VE p;
			int l = check_distant_mergable(x, y, ww, p);
			if(l < 0) continue;

			assert(l >= 2);
			assert(p.size() >= 2);
			assert(i2e[x] == p[0]);
			assert(i2e[y] == p[p.size() - 1]);

			bool c = false;
			for(int k = 1; k < p.size() - 1; k++)
			{
				int e = e2i[p[k]];
				if(ss.find(e) != ss.end()) c = true;
				if(c == true) break;
			}

			if(c == true) continue;

			vector<int> vv;
			split_merge_path(p, ww, vv);

			//// printf("connect (distant) edge pair (%d, %d)\n", x, y);

			md++;

			int n = vv.size() - 1;

			if(i2e[vv[0]] == null_edge) assert(vv[0] == x);
			if(i2e[vv[n]] == null_edge) assert(vv[n] == y);

			if(vv[0] == x) s.erase(s.begin() + i);
			else s[i] = vv[0];

			if(vv[n] == y) t.erase(t.begin() + j);
			else t[j] = vv[n];

			int f = resolve_equation(s, t, ma, md);
			if(f == 2) return 2;
			else return 1;
		}
	}

	return 0;
}

bool scallop2::check_adjacent_mergable(int ex, int ey, nested_graph &nt)
{
	assert(i2e[ex] != null_edge);
	assert(i2e[ey] != null_edge);

	int xs = i2e[ex]->source();
	int xt = i2e[ex]->target();
	int ys = i2e[ey]->source();
	int yt = i2e[ey]->target();

	if(xt == ys) return true;
	if(yt == xs) return true;

	bool b = false;
	vector<PI> xp, yp;
	if(gr.check_path(i2e[ex], i2e[ey])) b = nt.link(xs, xt, ys, yt, xp, yp);
	else if(gr.check_path(i2e[ey], i2e[ex])) b = nt.link(ys, yt, xs, xt, yp, xp);
	else return false;

	return b;
}

bool scallop2::check_adjacent_mergable(int ex, int ey, vector<PI> &p)
{
	assert(i2e[ex] != null_edge);
	assert(i2e[ey] != null_edge);

	int xs = i2e[ex]->source();
	int xt = i2e[ex]->target();
	int ys = i2e[ey]->source();
	int yt = i2e[ey]->target();

	if(xt == ys) return true;
	if(yt == xs) return true;

	vector<PI> xp, yp;
	bool b = false;

	nested_graph nt(gr);

	if(gr.check_path(i2e[ex], i2e[ey])) b = nt.link(xs, xt, ys, yt, xp, yp);
	else if(gr.check_path(i2e[ey], i2e[ex])) b = nt.link(ys, yt, xs, xt, yp, xp);
	else return false;
	
	if(b == false) return false;

	p = xp;
	p.insert(p.end(), yp.begin(), yp.end());

	return true;
}

int scallop2::check_distant_mergable(int x, int y, double w)
{
	VE p;
	return check_distant_mergable(x, y, w, p);
}

int scallop2::check_distant_mergable(int x, int y, double w, VE &p)
{
	p.clear();

	assert(i2e[x] != null_edge);
	assert(i2e[y] != null_edge);

	edge_descriptor xx = i2e[x];
	edge_descriptor yy = i2e[y];

	int xs = (xx)->source();
	int xt = (xx)->target();
	int ys = (yy)->source();
	int yt = (yy)->target();

	if(gr.check_path(yt, xs) == true)
	{
		int r = check_distant_mergable(y, x, w, p);
		reverse(p);
		return r;
	}

	if(gr.check_path(xt, ys) == false) return -1;

	int l = gr.compute_shortest_path_w(xt, ys, w, p);
	if(l < 0) return -1;

	p.insert(p.begin(), xx);
	p.insert(p.end(), yy);

	assert(p.size() == l + 2);
	return l + 2;
}

int scallop2::build_adjacent_edges(const vector<PI> &p)
{
	for(int i = 0; i < p.size(); i++)
	{
		int x = p[i].first;
		int y = p[i].second;
		if(y == -1)
		{
			int l = gr.compute_in_partner(x);
			int r = gr.compute_out_partner(x);
			gr.exchange(l, x, r);
		}
		else
		{
			gr.rotate(x, y);
		}
	}
	return 0;
}

bool scallop2::join_trivial_edges()
{
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		edge_descriptor e = *it1;
		int s = e->source();
		int t = e->target();

		if(s == 0) continue;
		if(t == gr.num_vertices() - 1) continue;
		if(gr.out_degree(s) != 1) continue;
		if(gr.in_degree(t) != 1) continue;

		bool b = join_trivial_edge(e);
		if(b == true) return true;
	}
	return false;
}

bool scallop2::join_trivial_edge(edge_descriptor &e)
{
	int s = e->source();
	int t = e->target();

	double ds = gr.get_vertex_stddev(s);
	double dt = gr.get_vertex_stddev(t);
	double as = gr.get_vertex_weight(s);
	double at = gr.get_vertex_weight(t);
	double dev = ds + dt;
	double ave = (as * ds + at * dt) / dev;
	gr.set_vertex_weight(s, ave);
	gr.set_vertex_stddev(s, dev);

	if(gr.locate(t) == 0)
	{
		double ew = 0;
		edge_iterator it1, it2;
		for(tie(it1, it2) = gr.out_edges(t); it1 != it2; it1++)
		{
			ew += gr.get_edge_weight(*it1);
		}
		gr.set_edge_weight(e, ew);

		decompose_trivial_vertex(t);
		return true;
	}
	else if(gr.locate(s) == 0)
	{
		double ew = 0;
		edge_iterator it1, it2;
		for(tie(it1, it2) = gr.in_edges(s); it1 != it2; it1++)
		{
			ew += gr.get_edge_weight(*it1);
		}
		gr.set_edge_weight(e, ew);
		decompose_trivial_vertex(s);
		return true;
	}

	return false;
}

bool scallop2::join_trivial_vertices()
{
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		bool b = join_trivial_vertex(i);
		if(b == true) return true;
	}
	return false;
}

bool scallop2::join_trivial_vertex(int i)
{
	if(gr.in_degree(i) != 1) return false;
	if(gr.out_degree(i) != 1) return false;

	edge_iterator it1, it2;
	tie(it1, it2) = gr.in_edges(i);
	edge_descriptor e1 = *it1;
	tie(it1, it2) = gr.out_edges(i);
	edge_descriptor e2 = *it1;

	if(e1->source() == 0) return false;
	if(e2->target() == gr.num_vertices() - 1) return false;

	double w = gr.get_vertex_weight(i);
	double d = 0;
	d += gr.get_vertex_stddev(i);
	d += gr.get_edge_stddev(e1);
	d += gr.get_edge_stddev(e2);

	gr.set_edge_weight(e1, w);
	gr.set_edge_weight(e2, w);
	gr.set_edge_stddev(e1, d / 2.0);
	gr.set_edge_stddev(e2, d / 2.0);

	decompose_trivial_vertex(i);
	return true;
}

bool scallop2::decompose_trivial_internal_vertices()
{
	bool flag = false;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.degree(i) == 0) continue;
		if(gr.in_degree(i) >= 2 && gr.out_degree(i) >= 2) continue;
		if(gr.locate(i) != 0) continue;
		decompose_trivial_vertex(i);
		flag = true;
	}
	return flag;
}

bool scallop2::decompose_trivial_external_vertices()
{
	bool flag = false;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.degree(i) == 0) continue;
		if(gr.in_degree(i) >= 2 && gr.out_degree(i) >= 2) continue;
		if(gr.locate(i) == 0) continue;
		decompose_trivial_vertex(i);
		flag = true;
	}
	return flag;
}

bool scallop2::decompose_trivial_vertex(int i)
{
	if(i == 0) return false;
	if(i == gr.num_vertices() - 1) return false;

	equation eqn(0);
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.in_edges(i); it1 != it2; it1++)
	{
		edge_descriptor e = *it1;
		int ei = e2i[e];
		eqn.s.push_back(ei);
	}
	for(tie(it1, it2) = gr.out_edges(i); it1 != it2; it1++)
	{
		edge_descriptor e = *it1;
		int ei = e2i[e];
		eqn.t.push_back(ei);
	}

	printf("decompose trivial vertex %d\n", i);
	eqn.print(i);

	resolve_equation(eqn);
	return true;
}

int scallop2::greedy_decompose()
{
	while(true)
	{
		VE v;
		vector<int> vv;
		double w = gr.compute_maximum_path_w(v);
		if(w <= 0.0) break;
		int e = split_merge_path(v, w, vv);
		collect_path(e);
	}
	return 0;
}

int scallop2::collect_existing_st_paths()
{
	for(int i = 0; i < i2e.size(); i++)
	{
		if(i2e[i] == null_edge) continue;
		if(i2e[i]->source() != 0) continue;
		if(i2e[i]->target() != gr.num_vertices() - 1) continue;
		collect_path(i);
	}
	return 0;
}

int scallop2::collect_path(int e)
{
	assert(i2e[e] != null_edge);
	assert(i2e[e]->source() == 0);
	assert(i2e[e]->target() == gr.num_vertices() - 1);

	assert(mev.find(i2e[e]) != mev.end());
	vector<int> v = mev[i2e[e]];
	sort(v.begin(), v.end());

	assert(v[0] == 0);
	assert(v[v.size() - 1] < gr.num_vertices() - 1);
	v.push_back(gr.num_vertices() - 1);

	path p;
	p.abd = gr.get_edge_weight(i2e[e]);
	p.v = v;
	paths.push_back(p);

	gr.remove_edge(i2e[e]);
	e2i.erase(i2e[e]);
	i2e[e] = null_edge;

	return 0;
}

int scallop2::print()
{
	int n = 0;
	for(int i = 0; i < gr.num_vertices(); i++) 
	{
		if(gr.degree(i) >= 1) n++;
	}

	int p1 = gr.compute_num_paths();
	int p2 = gr.compute_decomp_paths();
	printf("statistics: %lu edges, %d vertices, total %d paths, %d required\n", gr.num_edges(), n, p1, p2);
	printf("finish round %d\n\n", round);

	if(output_tex_files == true)
	{
		draw_splice_graph(name + "." + tostring(round) + ".tex");
		nested_graph nt(gr);
		nt.draw(name + "." + tostring(round) + ".nt.tex");
	}

	round++;

	return 0;
}

int scallop2::draw_splice_graph(const string &file) 
{
	MIS mis;
	char buf[10240];
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		double w = gr.get_vertex_weight(i);
		double d = gr.get_vertex_stddev(i);
		//string s = gr.get_vertex_string(i);
		//sprintf(buf, "%d:%.0lf:%s", i, w, s.c_str());
		sprintf(buf, "%d:%.1lf:%.0lf", i, w, d);
		mis.insert(PIS(i, buf));
	}

	MES mes;
	for(int i = 0; i < i2e.size(); i++)
	{
		if(i2e[i] == null_edge) continue;
		double w = gr.get_edge_weight(i2e[i]);
		double d = gr.get_edge_stddev(i2e[i]);
		sprintf(buf, "%d:%.1lf", i, w);
		mes.insert(PES(i2e[i], buf));
	}
	gr.draw(file, mis, mes, 4.0);
	return 0;
}

