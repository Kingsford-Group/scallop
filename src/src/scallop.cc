#include "scallop.h"
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

scallop::scallop()
{}

scallop::scallop(const string &s, splice_graph &g)
	: name(s), gr(g)
{
	round = 0;
	assert(gr.check_fully_connected() == true);
}

scallop::~scallop()
{}

int scallop::clear()
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

int scallop::save(scallop &sc)
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

int scallop::load(scallop &sc)
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

int scallop::assemble()
{
	int c = classify();
	if(c == TRIVIAL) return 0;

	if(algo == "basic") return assemble0();
	if(algo == "core") return assemble1();
	if(algo == "full") return assemble2();
	if(algo == "greedy") return greedy();
	return 0;
}

int scallop::classify()
{
	assert(gr.num_vertices() >= 2);
	if(gr.num_vertices() == 2) return TRIVIAL;

	string s;	
	int p0 = gr.compute_num_paths();
	int p1 = gr.num_edges() - gr.num_vertices() + 2;
	assert(p0 >= p1);

	bool b = (p0 == p1) ? true : false;

	printf("\nprocess %s %s\n", name.c_str(), b ? "TRIVIAL" : "NORMAL");

	if(p0 == p1) return TRIVIAL;
	else return NORMAL;
}

int scallop::assemble0()
{
	//if(output_tex_files == true) gr.draw(name + "." + tostring(round++) + ".tex");

	smoother sm(gr);
	sm.solve();

	//if(output_tex_files == true) gr.draw(name + "." + tostring(round++) + ".tex");

	//gr.round_weights();
	//remove_empty_edges();

	//if(output_tex_files == true) gr.draw(name + "." + tostring(round++) + ".tex");

	init_super_edges();
	reconstruct_splice_graph();
	gr.get_edge_indices(i2e, e2i);

	//print();
	//collect_existing_st_paths();
	//printf("%s scallop0 solution %lu paths\n", name.c_str(), paths.size());

	return 0;
}

int scallop::assemble1()
{
	assemble0();

	int f = iterate();

	collect_existing_st_paths();
	//print();

	printf("%s core solution %lu paths, iteration = %d\n", name.c_str(), paths.size(), f);

	return 0;
}

int scallop::assemble2()
{
	assemble0();

	int f = iterate();

	greedy_decompose();
	assert(gr.num_edges() == 0);

	//print();

	printf("%s full solution %lu paths, iteration = %d\n", name.c_str(), paths.size(), f);

	return 0;
}

int scallop::greedy()
{
	assemble0();

	greedy_decompose();
	assert(gr.num_edges() == 0);

	//print();

	printf("%s greedy solution %lu paths\n", name.c_str(), paths.size());

	return 0;
}

int scallop::iterate()
{
	print();
	while(true)
	{
		bool t = false;
		while(true)
		{
			bool b = decompose_trivial_vertices();
			if(b == true) t = true;
			if(b == false) break;
		}
		if(t == true) print();

		int f = decompose_with_equations();

		print();
		
		if(f == 0) continue;
		else return f;
	}
}

int scallop::decompose_with_equations()
{
	vector<equation> eqns;
	identify_equation1(eqns);

	if(eqns.size() == 0) return -2;

	sort(eqns.begin(), eqns.end(), equation_cmp1);

	for(int k = 0; k < eqns.size(); k++)
	{
		scallop sc;
		save(sc);
		decompose_single_equation(eqns[k]);
		load(sc);
	}

	sort(eqns.begin(), eqns.end(), equation_cmp2);

	for(int i = 0; i < eqns.size(); i++) 
	{
		printf(" "); 
		eqns[i].print(i);
	}

	if(eqns[0].b == false) return -1;

	decompose_single_equation(eqns[0]);

	return 0;
}

int scallop::decompose_single_equation(equation &eqn)
{
	vector<int> subs = eqn.s;
	vector<int> subt = eqn.t;

	VE vx, vy;
	for(int i = 0; i < subs.size(); i++) vx.push_back(i2e[subs[i]]);
	for(int i = 0; i < subt.size(); i++) vy.push_back(i2e[subt[i]]);

	smoother sm(gr);
	sm.add_equation(vx, vy);
	sm.solve();

	assert(subs.size() == 1);
	assert(subt.size() >= 1);

	subs = split_edge(subs[0], subt);

	PI p = connect_pairs(subs, subt);

	eqn.d = p.second;

	if(p.first + p.second >= subs.size()) eqn.b = true;
	else eqn.b = false;

	return 0;
}

int scallop::init_super_edges()
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

int scallop::reconstruct_splice_graph()
{
	while(true)
	{
		bool flag = false;
		for(int i = 0; i < gr.num_vertices(); i++)
		{
			bool b = init_trivial_vertex(i);
			if(b == true) flag = true;
		}
		if(flag == false) break;
	}
	return 0;
}

bool scallop::init_trivial_vertex(int x)
{
	int id = gr.in_degree(x);
	int od = gr.out_degree(x);

	if(id <= 0 || od <= 0) return false;
	if(id >= 2 && od >= 2) return false;
	//if(id <= 1 && od <= 1) return false;

	edge_iterator it1, it2;
	edge_iterator ot1, ot2;
	for(tie(it1, it2) = gr.in_edges(x); it1 != it2; it1++)
	{

		for(tie(ot1, ot2) = gr.out_edges(x); ot1 != ot2; ot1++)
		{
			int s = (*it1)->source();
			int t = (*ot1)->target();

			double w1 = gr.get_edge_weight(*it1);
			double a1 = gr.get_edge_stddev(*it1);
			double w2 = gr.get_edge_weight(*ot1);
			double a2 = gr.get_edge_stddev(*ot1);

			double w = w1 < w2 ? w1 : w2;
			double a = w1 < w2 ? a1 : a2;

			edge_descriptor p = gr.add_edge(s, t);
			gr.set_edge_weight(p, w);
			gr.set_edge_stddev(p, a);

			assert(mev.find(*it1) != mev.end());
			assert(mev.find(*ot1) != mev.end());

			vector<int> v1 = mev[*it1];
			vector<int> v2 = mev[*ot1];
			v1.insert(v1.end(), v2.begin(), v2.end());

			if(mev.find(p) != mev.end()) mev[p] = v1;
			else mev.insert(PEV(p, v1));
		}
	}
	gr.clear_vertex(x);
	return true;
}

int scallop::check_distant_mergable(int x, int y)
{
	VE p;
	return check_distant_mergable(x, y, p);
}

int scallop::check_distant_mergable(int x, int y, VE &p)
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

	double wx = gr.get_edge_weight(xx);
	double wy = gr.get_edge_weight(yy);
	assert(fabs(wx - wy) <= SMIN);

	if(gr.check_path(yt, xs) == true) return check_distant_mergable(y, x, p);
	if(gr.check_path(xt, ys) == false) return -1;

	int l = gr.compute_shortest_path_w(xt, ys, wx, p);
	if(l < 0) return -1;

	p.insert(p.begin(), xx);
	p.insert(p.end(), yy);

	return l + 2;
}

int scallop::connect_path(const VE &p, double wx)
{
	vector<int> v;
	for(int i = 0; i < p.size(); i++)
	{
		assert(p[i] != null_edge);
		assert(e2i.find(p[i]) != e2i.end());
		v.push_back(e2i[p[i]]);
	}
	return connect_path(v, wx);
}

int scallop::connect_path(const vector<int> &p, double ww)
{
	if(p.size() == 0) return -1;
	if(p.size() == 1) return p[0];

	int pret = i2e[p[0]]->target();
	for(int i = 1; i < p.size(); i++)
	{
		int s = i2e[p[i]]->source();
		int t = i2e[p[i]]->target();
		assert(s == pret);
		pret = t;
	}

	int pree = split_edge(p[0], ww);
	for(int i = 1; i < p.size(); i++)
	{
		int e = split_edge(p[i], pree);
		int ee = connect_adjacent_equal_edges(pree, e);
		pree = ee;
	}
	return pree;
}

int scallop::connect_adjacent_equal_edges(int x, int y)
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
	if(yt == xs) return connect_adjacent_equal_edges(y, x);
	
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
	gr.set_edge_stddev(p, wx1);

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


int scallop::split_edge(int ei, double w)
{
	assert(i2e[ei] != null_edge);
	edge_descriptor ee = i2e[ei];

	double ww = gr.get_edge_weight(ee);
	double dd = gr.get_edge_stddev(ee);

	if(fabs(ww - w) <= SMIN) return ei;

	assert(ww >= w + SMIN);

	int s = ee->source();
	int t = ee->target();

	edge_descriptor p1 = gr.add_edge(s, t);
	edge_descriptor p2 = gr.add_edge(s, t);

	gr.set_edge_weight(p1, w);
	gr.set_edge_weight(p2, ww - w);
	gr.set_edge_stddev(p1, dd);
	gr.set_edge_stddev(p2, dd);

	if(mev.find(p1) != mev.end()) mev[p1] = mev[ee];
	else mev.insert(PEV(p1, mev[ee]));

	if(mev.find(p2) != mev.end()) mev[p2] = mev[ee];
	else mev.insert(PEV(p2, mev[ee]));

	int n = i2e.size();
	i2e.push_back(p1);
	i2e.push_back(p2);
	e2i.insert(PEI(p1, n));
	e2i.insert(PEI(p2, n + 1));

	gr.remove_edge(ee);
	e2i.erase(ee);
	i2e[ei] = null_edge;

	return n;
}

int scallop::split_edge(int exi, int eyi)
{
	assert(i2e[exi] != null_edge);
	assert(i2e[eyi] != null_edge);
	edge_descriptor ex = i2e[exi];
	edge_descriptor ey = i2e[eyi];

	double wx = gr.get_edge_weight(ex);
	double wy = gr.get_edge_weight(ey);

	return split_edge(exi, wy);
}

vector<int> scallop::split_edge(int ei, const vector<int> &sub)
{
	vector<int> v;
	int x = ei;
	for(int i = 0; i < sub.size(); i++)
	{
		int y = split_edge(x, sub[i]);
		if(i == sub.size() - 1) assert(y == x);
		if(y == x) assert(i == sub.size() - 1);
		v.push_back(y);
		x = y + 1;
	}
	assert(v.size() == sub.size());
	return v;
}

bool scallop::verify_equation_nontrivial(const vector<int> &subs, const vector<int> &subt)
{
	set<edge_descriptor> fb;
	for(int i = 0; i < subs.size(); i++)
	{
		edge_descriptor &e = i2e[subs[i]];
		assert(e != null_edge);
		fb.insert(e);
	}

	vector<int> v, b;
	set<edge_descriptor> ss;

	gr.bfs(0, v, b, ss, fb);
	for(int i = 0; i < subt.size(); i++)
	{
		edge_descriptor &e = i2e[subt[i]];
		assert(e != null_edge);
		if(ss.find(e) == ss.end()) return false;
	}

	gr.bfs_reverse(gr.num_vertices() - 1, v, b, ss, fb);
	for(int i = 0; i < subt.size(); i++)
	{
		edge_descriptor &e = i2e[subt[i]];
		assert(e != null_edge);
		if(ss.find(e) == ss.end()) return false;
	}

	fb.clear();
	for(int i = 0; i < subt.size(); i++)
	{
		edge_descriptor &e = i2e[subt[i]];
		assert(e != null_edge);
		fb.insert(e);
	}

	gr.bfs(0, v, b, ss, fb);
	for(int i = 0; i < subs.size(); i++)
	{
		edge_descriptor &e = i2e[subs[i]];
		assert(e != null_edge);
		if(ss.find(e) == ss.end()) return false;
	}

	gr.bfs_reverse(gr.num_vertices() - 1, v, b, ss, fb);
	for(int i = 0; i < subs.size(); i++)
	{
		edge_descriptor &e = i2e[subs[i]];
		assert(e != null_edge);
		if(ss.find(e) == ss.end()) return false;
	}

	return true;
}

int scallop::identify_equation1(vector<equation> &eqns)
{
	eqns.clear();
	vector<PI> p;
	for(int i = 0; i < i2e.size(); i++)
	{
		if(i2e[i] == null_edge) continue;
		int w = (int)(gr.get_edge_weight(i2e[i]));
		if(w <= 0) continue;
		p.push_back(PI(w, i));
	}
	sort(p.begin(), p.end());

	for(int i = 0; i < p.size(); i++)
	{
		int e = p[i].second;
		assert(i2e[e] != null_edge);
		vector<int> s;
		s.push_back(e);
		vector<int> t;
		int err = identify_equation(s, t);

		double ratio = err * 1.0 / gr.get_edge_weight(i2e[e]);
		if(ratio > max_equation_error_ratio) continue;

		equation eqn(s, t, err);
		eqns.push_back(eqn);
	}
	return 0;
}

int scallop::identify_equation(const vector<int> &subs, vector<int> &subt)
{
	if(subs.size() == 0) return false;

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
		if(ss.find(e2i[*it]) != ss.end()) continue;
		int ww = (int)(gr.get_edge_weight(*it));
		if(ww <= 0) continue;
		if(ww * 1.0 > sw * 1.0 * (1.0 + max_equation_error_ratio)) continue;
		xi.push_back(PI(ww, e2i[*it]));
	}

	if(xi.size() == 0) return INT_MAX;

	sort(xi.begin(), xi.end());

	xi.push_back(PI(sw, e));

	vector<int> v;
	for(int i = 0; i < xi.size(); i++)
	{
		v.push_back(xi[i].first);
	}

	subsetsum sss(v);
	sss.solve();

	if(sss.subsets.size() == 0) return INT_MAX;

	for(int i = 0; i < sss.subsets.size(); i++)
	{
		vector<int> subset = sss.subsets[i];
		int opt = sss.opts[i];

		subt.clear();
		for(int j = 0; j < subset.size(); j++)
		{
			int k = subset[j];
			assert(xi[k].second != -1);
			subt.push_back(xi[k].second);
		}

		if(subt.size() >= 3) continue;
		if(verify_equation_nontrivial(subs, subt) == false) continue;

		int err = (int)fabs(opt - sw);
		return err;
	}
	
	return INT_MAX;
}

PI scallop::connect_pairs(const vector<int> &vx, const vector<int> &vy)
{
	PI pi(0, 0);
	assert(vx.size() == vy.size());

	vector<bool> m;
	m.assign(vx.size(), false);

	int n = (int)(gr.num_edges()) * 2;

	while(true)
	{
		vector<PI> r;
		for(int i = 0; i < vx.size(); i++)
		{
			r.push_back(PI(n, i));
		}

		for(int i = 0; i < vx.size(); i++)
		{
			int x = vx[i];
			int y = vy[i];

			if(i2e[x] == null_edge) continue;
			if(i2e[y] == null_edge) continue;

			double wx = gr.get_edge_weight(i2e[x]);
			double wy = gr.get_edge_weight(i2e[y]);

			if(fabs(wx - wy) > SMIN) continue;

			if(check_adjacent_mergable( x, y) == true)
			{
				r[i].first = 0;
				break;
			}

			int l = check_distant_mergable(x, y);
			assert(l < n);
			if(l < 0) continue;

			r[i].first = l;
			break;
		}

		sort(r.begin(), r.end());

		if(r[0].first == n) break;

		int k = r[0].second;
		m[k] = true;

		int x = vx[k];
		int y = vy[k];
		double w = gr.get_edge_weight(i2e[x]);

		if(r[0].first == 0)
		{
			vector<PI> p;
			bool b = check_adjacent_mergable(x, y, p);
			assert(b == true);
			build_adjacent_equal_edges(p);
			connect_adjacent_equal_edges(x, y);
			pi.first++;
			printf(" connect (adjacent) edge pair (%d, %d)\n", x, y);
		}
		else
		{
			VE p;
			int l = check_distant_mergable(x, y, p);
			assert(l >= 2);
			connect_path(p, w);
			pi.second++;
			printf(" connect (distant) edge pair (%d, %d)\n", x, y);
		}
	}

	return pi;
}

bool scallop::check_adjacent_mergable(int ex, int ey)
{
	vector<PI> p;
	return check_adjacent_mergable(ex, ey, p);
}

bool scallop::check_adjacent_mergable(int ex, int ey, vector<PI> &p)
{
	assert(i2e[ex] != null_edge);
	assert(i2e[ey] != null_edge);

	int xs = i2e[ex]->source();
	int xt = i2e[ex]->target();
	int ys = i2e[ey]->source();
	int yt = i2e[ey]->target();

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

int scallop::build_adjacent_equal_edges(const vector<PI> &p)
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

bool scallop::decompose_trivial_vertices()
{
	bool flag = false;
	edge_iterator it1, it2;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.degree(i) == 0) continue;
		if(gr.in_degree(i) == 1)
		{
			printf(" decompose trivial vertex %d\n", i);

			tie(it1, it2) = gr.in_edges(i);
			int ei = e2i[*it1];
			vector<int> sub;
			for(tie(it1, it2) = gr.out_edges(i); it1 != it2; it1++)
			{
				int e = e2i[*it1];
				sub.push_back(e);
			}

			vector<int> v = split_edge(ei, sub);
			assert(v.size() == sub.size());
			for(int k = 0; k < v.size(); k++)
			{
				assert(i2e[v[k]] != null_edge);
				connect_adjacent_equal_edges(v[k], sub[k]);
			}
			flag = true;
		}
		else if(gr.out_degree(i) == 1)
		{
			printf(" decompose trivial vertex %d\n", i);

			tie(it1, it2) = gr.out_edges(i);
			int ei = e2i[*it1];
			vector<int> sub;
			for(tie(it1, it2) = gr.in_edges(i); it1 != it2; it1++)
			{
				int e = e2i[*it1];
				sub.push_back(e);
			}

			vector<int> v = split_edge(ei, sub);
			assert(v.size() == sub.size());
			for(int k = 0; k < v.size(); k++)
			{
				assert(i2e[v[k]] != null_edge);
				connect_adjacent_equal_edges(sub[k], v[k]);
			}
			flag = true;
		}
	}
	return flag;
}

int scallop::greedy_decompose()
{
	while(true)
	{
		VE v;
		double w = gr.compute_maximum_path_w(v);
		if(w <= 0.0) break;
		int e = connect_path(v, w);
		collect_path(e);
	}
	return 0;
}

int scallop::collect_existing_st_paths()
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

int scallop::collect_path(int e)
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

int scallop::remove_empty_edges()
{
	for(int i = 0; i < i2e.size(); i++)
	{
		if(i2e[i] == null_edge) continue;
		double w = gr.get_edge_weight(i2e[i]);
		if(w >= 1) continue;
		assert(w <= 0);
		e2i.erase(i2e[i]);
		gr.remove_edge(i2e[i]);
		i2e[i] = null_edge;
	}
	return 0;
}

int scallop::print()
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

int scallop::draw_splice_graph(const string &file) 
{
	MIS mis;
	char buf[10240];
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		double w = gr.get_vertex_weight(i);
		sprintf(buf, "%d:%.0lf", i, w);
		mis.insert(PIS(i, buf));
	}

	MES mes;
	for(int i = 0; i < i2e.size(); i++)
	{
		if(i2e[i] == null_edge) continue;
		double w = gr.get_edge_weight(i2e[i]);
		sprintf(buf, "%d:%.0lf", i, w);
		//sprintf(buf, "%d", i);
		mes.insert(PES(i2e[i], buf));
	}
	gr.draw(file, mis, mes, 4.5);
	return 0;
}

