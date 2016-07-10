#include "scallop.h"
#include "subsetsum.h"
#include "subsetsum1.h"
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
	//assert(gr.check_fully_connected() == true);
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

int scallop::assemble0()
{
	if(output_tex_files == true) gr.draw(name + "." + tostring(round++) + ".tex");

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
		bool b = false;

		b = decompose_trivial_vertices();
		if(b == true) print();
		if(b == true) continue;

		b = decompose_with_equations();
		if(b == true) print();
		if(b == true) continue;

		break;
	}
	return 0;
}

bool scallop::decompose_with_equations()
{
	vector<equation> eqns;

	bool b = false;
	if(b == false) b = identify_equations1(eqns);
	if(b == false) b = identify_equations2(eqns);
	if(b == false) b = identify_equations3(eqns);

	if(eqns.size() == 0) return false;

	printf("equations = %lu\n", eqns.size());

	//resolve_vertex_with_equations(eqns);

	sort(eqns.begin(), eqns.end(), equation_cmp2);

	for(int i = 0; i < eqns.size(); i++) 
	{
		printf(" "); 
		eqns[i].print(i);
	}

	if(eqns[0].f == 3) return true;

	smooth_with_equation(eqns[0]);
	resolve_equation(eqns[0]);

	if(eqns[0].f >= 2) return true;
	else return false;
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

int scallop::split_merge_path(const VE &p, double wx, vector<int> &vv)
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

int scallop::split_merge_path(const vector<int> &p, double ww, vector<int> &vv)
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

int scallop::merge_adjacent_equal_edges(int x, int y)
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

int scallop::merge_adjacent_edges(int x, int y)
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

	edge_descriptor p2 = gr.add_edge(s, t);

	gr.set_edge_weight(ee, w);
	gr.set_edge_weight(p2, ww - w);
	gr.set_edge_stddev(p2, dd);

	if(mev.find(p2) != mev.end()) mev[p2] = mev[ee];
	else mev.insert(PEV(p2, mev[ee]));

	int n = i2e.size();
	i2e.push_back(p2);
	e2i.insert(PEI(p2, n));

	return n;
}

bool scallop::verify_equation_mergable(equation &eqn)
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

bool scallop::verify_equation_nontrivial(equation &eqn)
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
	if(b1 == false && b2 == false) return false;

	b1 = gr.bfs(vty, n, fbs);
	b2 = gr.bfs_reverse(vsx, 0, fbt);
	if(b1 == false && b2 == false) return false;

	return true;
}

bool scallop::identify_equations1(vector<equation> &eqns)
{
	typedef map<int, vector<int> > MIV;
	typedef pair<int, vector<int> > PIV;

	MIV m;
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		int w = (int)(gr.get_edge_weight(*it1));
		int e = e2i[*it1];
		if(m.find(w) != m.end()) 
		{
			m[w].push_back(e);
		}
		else
		{
			vector<int> v;
			v.push_back(e);
			m.insert(PIV(w, v));
		}
	}

	nested_graph nt(gr);
	for(MIV::iterator it = m.begin(); it != m.end(); it++)
	{
		vector<int> &v = it->second;
		if(v.size() <= 1) continue;
		for(int i = 0; i < v.size(); i++)
		{
			for(int j = i + 1; j < v.size(); j++)
			{
				bool b = check_adjacent_mergable(v[i], v[j], nt);
				int l = check_distant_mergable(v[i], v[j], it->first);

				vector<int> s;
				vector<int> t;
				s.push_back(v[i]);
				t.push_back(v[j]);

				equation eqn(s, t, 0);

				if(b == true || l >= 0) eqn.f = 2;
				else eqn.f = 1;

				if(b == true)
				{
					eqn.a = 1;
					eqn.d = 0;
				}
				else
				{
					eqn.a = 0;
					eqn.d = 1;
				}

				eqn.w = it->first;
				eqns.push_back(eqn);

				if(b == true) return true;
			}
		}
	}
	return false;
}

bool scallop::identify_equations2(vector<equation> &eqns)
{
	vector<PI> p;
	for(int i = 0; i < i2e.size(); i++)
	{
		edge_descriptor &e = i2e[i];
		if(e == null_edge) continue;
		if(e->source() == 0 && e->target() == gr.num_vertices() - 1) continue;
		int w = (int)(gr.get_edge_weight(e));
		if(w <= 0) continue;
		p.push_back(PI(w, i));
	}

	if(p.size() == 0) return false;

	subsetsum1 sss1(p, p);
	sss1.solve();

	for(int i = 0; i < sss1.eqns.size(); i++)
	{
		equation &eqn = sss1.eqns[i];
		assert(eqn.s.size() == 1);
		if(eqn.t.size() == 1) continue;
		if(verify_equation_mergable(eqn) == false) continue;
		if(verify_equation_nontrivial(eqn) == false) continue;

		scallop sc;
		save(sc);
		smooth_with_equation(eqn);
		resolve_equation(eqn);
		load(sc);

		eqns.push_back(eqn);
		if(eqn.f == 2 && eqn.d == 0) return true;
	}
	return false;
}

bool scallop::identify_equations3(vector<equation> &eqns)
{
	vector<PI> p;
	for(int i = 0; i < i2e.size(); i++)
	{
		edge_descriptor &e = i2e[i];
		if(e == null_edge) continue;
		if(e->source() == 0 && e->target() == gr.num_vertices() - 1) continue;
		int w = (int)(gr.get_edge_weight(e));
		if(w <= 0) continue;
		p.push_back(PI(w, i));
	}

	if(p.size() == 0) return false;

	vector<PI> ppi;
	vector<PI> ppw;
	for(int i = 0; i < p.size(); i++)
	{
		for(int j = i + 1; j < p.size(); j++)
		{
			int w = p[i].first + p[j].first;
			ppw.push_back(PI(w, ppw.size()));
			ppi.push_back(PI(p[i].second, p[j].second));
		}
	}

	subsetsum1 sss1(p, ppw);
	sss1.solve();

	for(int i = 0; i < sss1.eqns.size(); i++)
	{
		equation &eqn = sss1.eqns[i];
		assert(eqn.s.size() == 1);
		if(eqn.t.size() == 1) continue;

		int si = eqn.s[0];
		eqn.s.clear();
		eqn.s.push_back(ppi[si].first);
		eqn.s.push_back(ppi[si].second);

		set<int> ss(eqn.t.begin(), eqn.t.end());
		if(ss.find(eqn.s[0]) != ss.end()) continue;
		if(ss.find(eqn.s[1]) != ss.end()) continue;

		if(verify_equation_mergable(eqn) == false) continue;
		if(verify_equation_nontrivial(eqn) == false) continue;

		scallop sc;
		save(sc);
		smooth_with_equation(eqn);
		resolve_equation(eqn);
		load(sc);

		eqns.push_back(eqn);
		if(eqn.f == 2 && eqn.d == 0) return true;
	}

	/*
	nested_graph nt(gr);
	for(int i = 0; i < eqns.size(); i++)
	{
		equation & eqn = eqns[i];
		for(int j = 0; j < eqn.s.size(); j++)
		{
			for(int k = 0; k < eqn.t.size(); k++)
			{
				bool b = check_adjacent_mergable(eqn.s[j], eqn.t[k], nt);
				if(b == false) eqn.d++;
				else eqn.a++;
			}
		}
	}
	*/
	return false;
}

bool scallop::resolve_vertex_with_equations(vector<equation> &eqns)
{
	for(int i = 0; i < eqns.size(); i++)
	{
		resolve_vertex_with_equation1(eqns[i]);
		resolve_vertex_with_equation2(eqns[i]);
		if(eqns[i].f == 3) return true;
	}
	return false;
}

bool scallop::resolve_vertex_with_equation1(equation &eqn)
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
		merge_adjacent_equal_edges(ei, e2i[e]);
		ei = ee;
	}

	eqn.f = 3;
	return true;
}

bool scallop::resolve_vertex_with_equation2(equation &eqn)
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
		merge_adjacent_equal_edges(ei, e2i[e]);
		ei = ee;
	}

	eqn.f = 3;
	return true;
}

int scallop::smooth_with_equation(equation &eqn)
{
	return 0;		// TODO

	VE vx, vy;
	for(int i = 0; i < eqn.s.size(); i++) vx.push_back(i2e[eqn.s[i]]);
	for(int i = 0; i < eqn.t.size(); i++) vy.push_back(i2e[eqn.t[i]]);

	smoother sm(gr);
	sm.add_equation(vx, vy);
	sm.solve();
	
	return 0;
}

int scallop::resolve_equation(equation &eqn)
{
	vector<int> s = eqn.s;
	vector<int> t = eqn.t;
	eqn.a = eqn.d = 0;
	eqn.f = resolve_equation(s, t, eqn.a, eqn.d);
	return 0;
}

int scallop::resolve_equation(vector<int> &s, vector<int> &t, int &ma, int &md)
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

bool scallop::check_adjacent_mergable(int ex, int ey, nested_graph &nt)
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

bool scallop::check_adjacent_mergable(int ex, int ey, vector<PI> &p)
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

int scallop::check_distant_mergable(int x, int y, double w)
{
	VE p;
	return check_distant_mergable(x, y, w, p);
}

int scallop::check_distant_mergable(int x, int y, double w, VE &p)
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

int scallop::build_adjacent_edges(const vector<PI> &p)
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
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.degree(i) == 0) continue;
		if(gr.in_degree(i) >= 2 && gr.out_degree(i) >= 2) continue;

		printf("decompose trivial vertex %d\n", i);

		equation eqn(0);
		edge_iterator it1, it2;
		for(tie(it1, it2) = gr.in_edges(i); it1 != it2; it1++)
		{
			int e = e2i[*it1];
			eqn.s.push_back(e);
		}
		for(tie(it1, it2) = gr.out_edges(i); it1 != it2; it1++)
		{
			int e = e2i[*it1];
			eqn.t.push_back(e);
		}

		printf(" ");
		eqn.print(i);

		resolve_equation(eqn);

		flag = true;
	}
	return flag;
}

int scallop::greedy_decompose()
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

