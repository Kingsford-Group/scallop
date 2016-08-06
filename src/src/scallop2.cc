#include "scallop2.h"
#include "subsetsum.h"
#include "config.h"
#include "smoother.h"
#include "nested_graph.h"
#include "undirected_graph.h"

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
	if(output_tex_files == true) gr.draw(name + "." + tostring(round++) + ".tex");
	gr.get_edge_indices(i2e, e2i);
	init_super_edges();
	//assert(gr.check_fully_connected() == true);
}

scallop2::scallop2(const string &s, splice_graph &g, const vector<hyper_edge> &vhe)
	: name(s), gr(g)
{
	round = 0;
	if(output_tex_files == true) gr.draw(name + "." + tostring(round++) + ".tex");
	gr.get_edge_indices(i2e, e2i);
	init_super_edges();
	init_gateways(vhe);
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

int scallop2::assemble()
{
	classify();
	if(algo == "basic") return assemble0();
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
	compute_shortest_source_distances();
	compute_shortest_target_distances();
	assign_reliability();

	print();

	infer_equivalent_classes();

	while(true)
	{
		bool f = false;
		bool b = false;
		while(true)
		{
			b = infer_trivial_edges();
			if(b == true) f = true;
			if(b == true) continue;

			b = infer_vertices();
			if(b == true) f = true;
			if(b == true) continue;

			break;
		}

		if(f == true) print();

		f = false;
		while(true)
		{
			b = join_trivial_edges();
			if(b == true) f = true;
			if(b == true) continue;

			b = join_trivial_vertices();
			if(b == true) f = true;
			if(b == true) continue;

			break;
		}

		if(f == true) print();
		if(f == true) continue;

		f = false;
		while(true)
		{
			b = infer_edges();
			if(b == true) f = true;
			if(b == true) continue;
			break;
		}

		if(f == true) print();
		if(f == true) continue;

		b = decompose_vertices();
		if(b == true) print();
		if(b == true) continue;

		break;
	}

	smooth_splice_graph();
	print();

	return 0;
}

int scallop2::assemble1()
{
	assemble0();

	iterate(false);

	collect_existing_st_paths();

	greedy_decompose(1);

	//iterate(false);
	//collect_existing_st_paths();

	printf("%s core solution %lu paths\n", name.c_str(), paths.size());

	return 0;
}

int scallop2::assemble2()
{
	assemble0();

	iterate(true);

	//greedy_decompose(-1);
	assert(gr.num_edges() == 0);

	printf("%s full solution %lu paths\n", name.c_str(), paths.size());

	return 0;
}

int scallop2::greedy()
{
	assemble0();

	greedy_decompose(-1);
	assert(gr.num_edges() == 0);

	printf("%s greedy solution %lu paths\n", name.c_str(), paths.size());

	return 0;
}

int scallop2::iterate(bool greedy)
{
	while(true)
	{
		bool b = false;

		b = remove_false_boundary_edges();
		if(b == true) smooth_splice_graph();
		if(b == true) print();
		if(b == true) continue;

		b = decompose_trivial_vertices();
		if(b == true) print();
		if(b == true) continue;

		b = decompose_with_equations(0);
		if(b == true) print();
		if(b == true) continue;

		b = decompose_with_equations(1);
		if(b == true) print();
		if(b == true) continue;

		if(greedy == false) break;

		collect_existing_st_paths();
		greedy_decompose(1);

		if(gr.num_edges() == 0) break;
	}

	return 0;
}


bool scallop2::decompose_with_equations(int level)
{
	vector<equation> eqns;

	//gr.round_weights();

	bool b = false;
	if(level == 0) b = identify_equations0(eqns);
	if(level == 1) b = identify_equations1(eqns);
	if(level == 2) b = identify_equations2(eqns);

	if(eqns.size() == 0) return false;

	sort(eqns.begin(), eqns.end(), equation_cmp1);

	printf("candidate equations: %lu\n", eqns.size());

	for(int i = 0; i < eqns.size(); i++)
	{
		eqns[i].print(i);
	}

	equation eqn(0);
	for(int i = 0; i < eqns.size(); i++)
	{
		eqn = eqns[i];

		if(eqn.f == 2 && eqn.d == 0) break;

		if(verify_equation_mergable(eqn) == false) continue;
		if(verify_equation_nontrivial(eqn) == false) continue;

		scallop2 sc;
		save(sc);

		bool b = smooth_splice_graph(eqn);
		if(b == true) resolve_equation(eqn);
		load(sc);

		if(eqn.f == 2 && eqn.d == 0) break;
	}

	if(eqn.f != 2 || eqn.d != 0) return false;

	printf("smooth with equation\n");
	eqn.print(99);
	smooth_splice_graph(eqn);
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

int scallop2::init_gateways(const vector<hyper_edge> &vhe)
{
	gateways.clear();
	gateways.resize(gr.num_vertices());

	for(int i = 0; i < vhe.size(); i++)
	{
		const hyper_edge &he = vhe[i];
		const vector<int> &vv = he.v;
		if(vv.size() <= 1) continue;

		VE ve;
		for(int k = 0; k < vv.size() - 1; k++)
		{
			PEB p = gr.edge(vv[k], vv[k + 1]);
			if(p.second == false) ve.push_back(null_edge);
			else ve.push_back(p.first);
		}

		for(int k = 0; k < ve.size() - 1; k++)
		{
			if(ve[k] == null_edge || ve[k + 1] == null_edge) continue;

			int e1 = e2i[ve[k]];
			int e2 = e2i[ve[k + 1]];
			int c = he.count;
			int x = vv[k + 1];
			gateways[x].add_route(PI(e1, e2), c);
		}
	}
	return 0;
}


bool scallop2::decompose_vertices()
{
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		if(vi.infer == false && vi.reliability <= infer_min_reliability) continue;
		if(gr.in_degree(i) <= 1) continue;
		if(gr.out_degree(i) <= 1) continue;
		bool b = decompose_vertex(i);
		if(b == true) return true;
	}
	return false;
}

bool scallop2::decompose_vertex(int x)
{
	gateway &gw = gateways[x];

	double c = gw.total_counts();
	if(c < min_hyper_edges_count) return false;

	undirected_graph ug;
	MI e2u;
	vector<int> u2e;
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.in_edges(x); it1 != it2; it1++)
	{
		if(gr.get_edge_info(*it1).length == 0 && (*it1)->source() == 0) return false;
		ug.add_vertex();
		int e = e2i[*it1];
		e2u.insert(PI(e, e2u.size()));
		u2e.push_back(e);
	}
	for(tie(it1, it2) = gr.out_edges(x); it1 != it2; it1++)
	{
		if(gr.get_edge_info(*it1).length == 0 && (*it1)->target() == gr.num_vertices() - 1) return false;
		ug.add_vertex();
		int e = e2i[*it1];
		e2u.insert(PI(e, e2u.size()));
		u2e.push_back(e);
	}

	for(int i = 0; i < gw.routes.size(); i++)
	{
		int e1 = gw.routes[i].first;
		int e2 = gw.routes[i].second;

		assert(e2u.find(e1) != e2u.end());
		assert(e2u.find(e2) != e2u.end());

		int s = e2u[e1];
		int t = e2u[e2];
		assert(s >= 0 && s < gr.in_degree(x));
		assert(t >= gr.in_degree(x) && t < gr.degree(x));

		ug.add_edge(s, t);
	}

	vector< set<int> > vv = ug.compute_connected_components();

	int k = -1;
	equation eqn0;
	double r0 = 1.0;
	for(int i = 0; i < vv.size(); i++)
	{
		int kk = -1;
		equation eqn;
		for(set<int>::iterator it = vv[i].begin(); it != vv[i].end(); it++)
		{
			int v = (*it);
			int e = u2e[v];
			if(ug.degree(v) == 1) kk = v;
			if(v < gr.in_degree(x)) eqn.s.push_back(e);
			else eqn.t.push_back(e);
		}

		if(eqn.s.size() <= 0 || eqn.t.size() <= 0) continue;
		if(kk == -1) continue;

		double r = compute_equation_error_ratio(eqn);
		if(r <= -1 || r >= max_gateway_error_ratio) continue;
		if(r > r0) continue;

		eqn0 = eqn;
		r0 = r;
		k = kk;
	}

	if(k == -1) return false;

	bool b = smooth_vertex(x, eqn0);
	if(b == false) return false;

	printf("smooth vertex %d with equation\n", x);
	eqn0.print(99);
	gw.print(x);

	assert(ug.degree(k) == 1);
	tie(it1, it2) = ug.out_edges(k);
	int t = (*it1)->target();
	assert(k != t);
	int xx = u2e[k];
	int yy = u2e[t];

	printf("decompose vertex %d\n", x);
	printf("merge adjacent edges (%d, %d)\n", xx, yy);

	int ee = merge_adjacent_edges(xx, yy);
	edge_info ei = gr.get_edge_info(i2e[ee]);
	ei.infer = true;
	gr.set_edge_info(i2e[ee], ei);

	gw.print(x);

	return true;
}

int scallop2::split_merge_path(const vector<int> &p, double ww)
{
	vector<int> vv;
	return split_merge_path(p, ww, vv);
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
	assert(fabs(wx0 - wy0) <= SMIN);

	int lx1 = gr.get_edge_info(xx).length;
	int ly1 = gr.get_edge_info(yy).length;
	int lxt = gr.get_vertex_info(xt).length;
	int lxy = lx1 + ly1 + lxt;

	gr.set_edge_weight(p, wx0);
	gr.set_edge_info(p, edge_info(lxy));

	double wv = gr.get_vertex_weight(xt);
	gr.set_vertex_weight(xt, wv - wx0);

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

	int xs = xx->source();
	int xt = xx->target();
	int ys = yy->source();
	int yt = yy->target();

	if(xt != ys) return merge_adjacent_edges(y, x);
	assert(xt == ys);

	double wx = gr.get_edge_weight(xx);
	double wy = gr.get_edge_weight(yy);
	double ww = (wx <= wy) ? wx : wy;

	int x1 = split_edge(x, ww);
	int y1 = split_edge(y, ww);

	gateways[xs].split_out_edge(x, x1, ww / wx);
	gateways[yt].split_in_edge(y, y1, ww / wy);

	int xy = merge_adjacent_equal_edges(x, y);

	gateways[xs].replace_out_edge(x, xy);
	gateways[yt].replace_in_edge(y, xy);
	gateways[xt].remove_in_edge(x);
	gateways[xt].remove_out_edge(y);

	return xy;
}

int scallop2::split_edge(int ei, double w)
{
	assert(i2e[ei] != null_edge);
	edge_descriptor ee = i2e[ei];

	double ww = gr.get_edge_weight(ee);

	if(fabs(ww - w) <= SMIN) return ei;
	assert(ww >= w + SMIN);

	int s = ee->source();
	int t = ee->target();

	edge_descriptor p2 = gr.add_edge(s, t);
	edge_info eif = gr.get_edge_info(ee);

	gr.set_edge_weight(ee, w);
	gr.set_edge_info(ee, eif);
	gr.set_edge_weight(p2, ww - w);
	gr.set_edge_info(p2, eif);

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
		int ww = (int)(gr.get_edge_weight(*it));
		if(ww <= 0) continue;
		if(ww > sw) continue;
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

		if(err * 1.0 > 1.0 * w * max_equation_error_ratio) continue;

		equation eqn(err);
		eqn.s = subs;
		eqn.t = subt;
		eqns.push_back(eqn);

		//return 0;		// only consider one equation
	}
	return 0;
}

int scallop2::identify_equations0(vector<equation> &eqns)
{
	typedef pair<double, int> PDI;
	vector<PDI> vv;
	for(int i = 0; i < i2e.size(); i++)
	{
		edge_descriptor &e = i2e[i];
		if(e == null_edge) continue;
		double w = gr.get_edge_weight(e);
		vv.push_back(PDI(w, i));
	}

	sort(vv.begin(), vv.end());

	nested_graph nt(gr);
	for(int i = vv.size() - 1; i >= 1; i--)
	{
		int i1 = vv[i].second;
		double w1 = vv[i].first;
		edge_descriptor &e1 = i2e[i1];
		int s1 = e1->source();
		int t1 = e1->target();

		SE ss, tt;
		gr.bfs(t1, tt);
		gr.bfs_reverse(s1, ss);

		for(int j = i - 1; j >= 0; j--)
		{
			int i2 = vv[j].second;
			double w2 = vv[j].first;
			edge_descriptor &e2 = i2e[i2];
			int s2 = e2->source();
			int t2 = e2->target();

			assert(w1 >= w2);

			if(ss.find(e2) != ss.end())
			{
				bool b = nt.link(s2, t2, s1, t1);
				if(b == false) continue;
			}
			else if(tt.find(e2) != tt.end())
			{
				bool b = nt.link(s1, t1, s2, t2);
				if(b == false) continue;
			}
			else continue;

			double err = fabs(w1 - w2);
			//printf("w1 = %.2lf, w2 = %.2lf, error = %.2lf, ratio = %.2lf\n", w1, w2, err, max_equation_error_ratio);
			if(err * 1.0 >= 1.0 * w1 * max_equation_error_ratio) continue;

			equation eqn(err);
			eqn.s.push_back(i1);
			eqn.t.push_back(i2);
			eqn.f = 2;
			eqn.d = 0;
			eqns.push_back(eqn);

			break;
		}
	}
	return 0;
}

int scallop2::identify_equations1(vector<equation> &eqns)
{
	for(int i = 0; i < i2e.size(); i++)
	{
		edge_descriptor &e = i2e[i];
		if(e == null_edge) continue;

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

bool scallop2::smooth_splice_graph()
{
	vector<equation> eqns;
	return smooth_splice_graph(eqns);
}

bool scallop2::smooth_splice_graph(const equation &eqn)
{
	vector<equation> eqns;
	eqns.push_back(eqn);
	return smooth_splice_graph(eqns);
}

bool scallop2::smooth_splice_graph(const vector<equation> &eqns)
{
	smoother sm(gr);
	for(int k = 0; k < eqns.size(); k++)
	{
		const equation &eqn = eqns[k];

		VE vx, vy;
		for(int i = 0; i < eqn.s.size(); i++) vx.push_back(i2e[eqn.s[i]]);
		for(int i = 0; i < eqn.t.size(); i++) vy.push_back(i2e[eqn.t[i]]);

		sm.add_equation(vx, vy);
	}

	int f = sm.smooth();

	if(f == 0) return true;
	else return false;
}

bool scallop2::smooth_vertex(int v)
{
	vector<equation> eqns;
	return smooth_vertex(v, eqns);
}

bool scallop2::smooth_vertex(int v, const equation &eqn)
{
	vector<equation> eqns;
	eqns.push_back(eqn);
	return smooth_vertex(v, eqns);
}

bool scallop2::smooth_vertex(int v, const vector<equation> &eqns)
{
	smoother sm(gr);
	for(int k = 0; k < eqns.size(); k++)
	{
		const equation &eqn = eqns[k];

		VE vx, vy;
		for(int i = 0; i < eqn.s.size(); i++) vx.push_back(i2e[eqn.s[i]]);
		for(int i = 0; i < eqn.t.size(); i++) vy.push_back(i2e[eqn.t[i]]);

		sm.add_equation(vx, vy);
	}

	int f = sm.smooth_vertex(v);

	if(f == 0) return true;
	else return false;
}

bool scallop2::smooth_vertex0(int v, const vector<equation> &eqns)
{
	if(gr.get_vertex_info(v).infer == false && gr.get_vertex_info(v).reliability <= infer_min_reliability) return false;
	
	map<int, bool> m1, m2;
	double sum1 = 0, sum2 = 0;
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.in_edges(v); it1 != it2; it1++)
	{
		int e = e2i[*it1];
		m1.insert(pair<int, bool>(e, false));
		sum1 += gr.get_edge_weight(*it1);
	}
	for(tie(it1, it2) = gr.out_edges(v); it1 != it2; it1++)
	{
		int e = e2i[*it1];
		m2.insert(pair<int, bool>(e, false));
		sum2 += gr.get_edge_weight(*it1);
	}

	if(sum1 <= SMIN) return false;
	if(sum2 <= SMIN) return false;

	vector<double> vv1, vv2;
	for(int i = 0; i < eqns.size(); i++)
	{
		const equation &eqn = eqns[i];
		double s1 = 0, s2 = 0;
		for(int k = 0; k < eqn.s.size(); k++)
		{
			int e = eqn.s[k];
			assert(m1.find(e) != m1.end());
			assert(m1[e] == false);
			m1[e] = true;
			s1 += gr.get_edge_weight(i2e[e]);
		}
		for(int k = 0; k < eqn.t.size(); k++)
		{
			int e = eqn.t[k];
			assert(m2.find(e) != m2.end());
			assert(m2[e] == false);
			m2[e] = true;
			s2 += gr.get_edge_weight(i2e[e]);
		}
		vv1.push_back(s1);
		vv2.push_back(s2);
	}

	vector<int> v1, v2;
	double s1 = 0, s2 = 0;
	for(map<int, bool>::iterator it = m1.begin(); it != m1.end(); it++)
	{
		if(it->second == true) continue;
		v1.push_back(it->first);
		s1 += gr.get_edge_weight(i2e[it->first]);
	}
	for(map<int, bool>::iterator it = m2.begin(); it != m2.end(); it++)
	{
		if(it->second == true) continue;
		v2.push_back(it->first);
		s2 += gr.get_edge_weight(i2e[it->first]);
	}

	double w = gr.get_vertex_weight(v);
	for(int i = 0; i < eqns.size(); i++)
	{
		double ww = (vv1[i] / sum1 + vv2[i] / sum2) * 0.5 * w;
		distribute_weights(ww, eqns[i].s);
		distribute_weights(ww, eqns[i].t);
	}

	double ww = (s1 / sum1 + s2 / sum2) * 0.5 * w;
	distribute_weights(ww, v1);
	distribute_weights(ww, v2);

	return true;
}

bool scallop2::smooth_trivial_equation(const equation &e)
{
	equation eqn = e;
	if(eqn.s.size() >= 2 && eqn.t.size() == 1)
	{
		vector<int> v = eqn.s;
		eqn.s = eqn.t;
		eqn.t = v;
		return smooth_trivial_equation(eqn);
	}

	if(eqn.s.size() != 1) return false;
	if(eqn.t.size() <= 0) return false;

	vector<double> cc;
	vector<double> ww;
	for(int i = 0; i < eqn.s.size(); i++)
	{
		edge_descriptor e = i2e[eqn.s[i]];
		double w = gr.get_edge_weight(e);
		double c = gr.get_edge_info(e).length + pseudo_length_count;
		ww.push_back(w);
		cc.push_back(c);
	}
	for(int i = 0; i < eqn.t.size(); i++)
	{
		edge_descriptor e = i2e[eqn.t[i]];
		double w = gr.get_edge_weight(e);
		double c = gr.get_edge_info(e).length + pseudo_length_count;
		ww.push_back(w);
		cc.push_back(c);
	}

	vector<double> ww1;

	double w0 = 0, c0 = 1.0;
	for(int i = 1; i < ww.size(); i++) c0 += (cc[0] / cc[i]);
	for(int i = 1; i < ww.size(); i++) w0 += (ww[i] + ww[0] * cc[0] / cc[i]);

	ww1.push_back(w0 / c0);

	for(int i = 1; i < ww.size(); i++)
	{
		double w = (cc[0] * ww[0] + cc[i] * ww[i] - cc[0] * ww1[0] ) / cc[i];
		ww1.push_back(w);
	}

	for(int i = 0; i < eqn.s.size(); i++)
	{
		edge_descriptor e = i2e[eqn.s[i]];
		gr.set_edge_weight(e, ww1[i]);
	}

	for(int i = 0; i < eqn.t.size(); i++)
	{
		edge_descriptor e = i2e[eqn.t[i]];
		gr.set_edge_weight(e, ww1[i]);
	}

	return 0;
}

int scallop2::resolve_equation(equation &eqn)
{
	vector<int> ve;
	return resolve_equation(eqn, ve);
}

int scallop2::resolve_equation(equation &eqn, vector<int> &ve)
{
	vector<int> s = eqn.s;
	vector<int> t = eqn.t;
	eqn.a = eqn.d = 0;
	eqn.f = resolve_equation(s, t, eqn.a, eqn.d, ve);
	return 0;
}

int scallop2::resolve_equation(vector<int> &s, vector<int> &t, int &ma, int &md, vector<int> &ve)
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
			int ee = split_merge_path(v, ww, vv);
			ve.push_back(ee);

			ma++;

			//// printf("merge (adjacent) edge pair (%d, %d)\n", x, y);

			if(i2e[vv[0]] == null_edge) assert(vv[0] == x);
			if(i2e[vv[1]] == null_edge) assert(vv[1] == y);

			if(vv[0] == x) s.erase(s.begin() + i);
			else s[i] = vv[0];

			if(vv[1] == y) t.erase(t.begin() + j);
			else t[j] = vv[1];

			int f = resolve_equation(s, t, ma, md, ve);
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
			int ee = split_merge_path(p, ww, vv);
			ve.push_back(ee);

			//// printf("connect (distant) edge pair (%d, %d)\n", x, y);

			md++;

			int n = vv.size() - 1;

			if(i2e[vv[0]] == null_edge) assert(vv[0] == x);
			if(i2e[vv[n]] == null_edge) assert(vv[n] == y);

			if(vv[0] == x) s.erase(s.begin() + i);
			else s[i] = vv[0];

			if(vv[n] == y) t.erase(t.begin() + j);
			else t[j] = vv[n];

			int f = resolve_equation(s, t, ma, md, ve);
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

	bool bs = gr.get_vertex_info(s).infer;
	bool bt = gr.get_vertex_info(t).infer;
	double r1 = gr.get_vertex_info(s).reliability;
	double r2 = gr.get_vertex_info(t).reliability;

	if(bs == false && r1 <= join_min_reliability) return false;
	if(bt == false && r2 <= join_min_reliability) return false;

	double ws = gr.get_vertex_weight(s);
	double wt = gr.get_vertex_weight(t);
	double we = gr.get_edge_weight(e);

	int ls = gr.get_vertex_info(s).length;
	int lt = gr.get_vertex_info(t).length;
	int le = gr.get_edge_info(e).length;

	int lsp = ls + pseudo_length_count;
	int ltp = lt + pseudo_length_count;
	int lep = le + pseudo_length_count;

	if(ws <= SMIN) lsp = 0;
	if(wt <= SMIN) ltp = 0;
	if(we <= SMIN) lep = 0;

	double ww = (lsp * ws + ltp * wt + lep * we) / (lsp + ltp + lep);

	if(gr.locate(t) == 0)
	{
		gr.set_vertex_weight(s, ww);
		vertex_info vi = gr.get_vertex_info(s);
		vi.length = ls + lt + le;
		gr.set_vertex_info(s, vi);

		double ew = 0;
		edge_iterator it1, it2;
		for(tie(it1, it2) = gr.out_edges(t); it1 != it2; it1++)
		{
			ew += gr.get_edge_weight(*it1);
		}
		gr.set_edge_weight(e, ew);
		edge_info ei = gr.get_edge_info(e);
		vi = gr.get_vertex_info(t);
		ei.length = 0;
		vi.length = 0;
		gr.set_edge_info(e, 0);
		gr.set_vertex_info(t, vi);

		printf("join trivial edge %d\n", e2i[e]);
		decompose_trivial_vertex(t);
		return true;
	}
	else if(gr.locate(s) == 0)
	{
		gr.set_vertex_weight(t, ww);
		vertex_info vi = gr.get_vertex_info(t);
		vi.length = ls + lt + le;
		gr.set_vertex_info(t, vi);

		double ew = 0;
		edge_iterator it1, it2;
		for(tie(it1, it2) = gr.in_edges(s); it1 != it2; it1++)
		{
			ew += gr.get_edge_weight(*it1);
		}
		gr.set_edge_weight(e, ew);
		edge_info ei = gr.get_edge_info(e);
		vi = gr.get_vertex_info(s);
		ei.length = 0;
		vi.length = 0;
		gr.set_edge_info(e, ei);
		gr.set_vertex_info(s, vi);

		printf("join trivial edge %d\n", e2i[e]);
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

	double r = gr.get_vertex_info(i).reliability;
	bool b = gr.get_vertex_info(i).infer;
	if(b == false && r <= join_min_reliability) return false;

	edge_iterator it1, it2;
	tie(it1, it2) = gr.in_edges(i);
	edge_descriptor e1 = *it1;
	tie(it1, it2) = gr.out_edges(i);
	edge_descriptor e2 = *it1;

	//if(e1->source() == 0) return false;
	//if(e2->target() == gr.num_vertices() - 1) return false;

	double wv = gr.get_vertex_weight(i);
	double w1 = gr.get_edge_weight(e1);
	double w2 = gr.get_edge_weight(e2);
	int lv = gr.get_vertex_info(i).length + pseudo_length_count;
	int l1 = gr.get_edge_info(e1).length + pseudo_length_count;
	int l2 = gr.get_edge_info(e2).length + pseudo_length_count;

	if(w1 <= SMIN) l1 = 0;
	if(w2 <= SMIN) l2 = 0;

	double w = (wv * lv + w1 * l1 + w2 * l2) / (lv + l1 + l2);

	gr.set_edge_weight(e1, w);
	gr.set_edge_weight(e2, w);

	printf("join trivial vertex %d\n", i);
	vector<int> ve;
	decompose_trivial_vertex(i, ve);
	assert(ve.size() == 1);
	edge_descriptor e = i2e[ve[0]];
	edge_info ei = gr.get_edge_info(e);
	ei.infer = gr.get_vertex_info(i).infer;
	gr.set_edge_info(e, ei);

	return true;
}

int scallop2::compute_shortest_source_distances()
{
	MED med = gr.get_edge_weights();
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int le = gr.get_edge_info(*it1).length;
		int lv = gr.get_vertex_info(s).length;
		gr.set_edge_weight(*it1, le + lv);
	}

	vector<double> d;
	gr.compute_closest_path(0, d);

	assert(d.size() == gr.num_vertices());

	for(int i = 0; i < d.size(); i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		vi.sdist = d[i];
		gr.set_vertex_info(i, vi);
	}

	gr.set_edge_weights(med);
	return 0;
}

int scallop2::compute_shortest_target_distances()
{
	MED med = gr.get_edge_weights();
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		int t = (*it1)->target();
		int le = gr.get_edge_info(*it1).length;
		int lv = gr.get_vertex_info(t).length;
		gr.set_edge_weight(*it1, le + lv);
	}

	vector<double> d;
	gr.compute_closest_path_reverse(gr.num_vertices() - 1, d);

	assert(d.size() == gr.num_vertices());

	for(int i = 0; i < d.size(); i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		vi.tdist = d[i];
		gr.set_vertex_info(i, vi);
	}

	gr.set_edge_weights(med);

	return 0;
}

int scallop2::assign_reliability()
{
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		int d = vi.sdist < vi.tdist ? vi.sdist : vi.tdist;
		double r = d * 1.0 / average_slope_length;
		if(r >= 1.0) r = 1.0;
		if(r <= 0.0) r = 0.0;
		vi.reliability = r;
		gr.set_vertex_info(i, vi);
	}
	return 0;
}

bool scallop2::infer_equivalent_classes()
{
	vector<bool> m;
	m.assign(gr.num_vertices(), false);

	bool flag = false;
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		if(m[i] == true) continue;

		vector<int> v;
		int x = i;
		while(x != -1)
		{
			assert(m[x] == false);
			m[x] = true;
			v.push_back(x);
			x = gr.compute_out_equivalent_vertex(x);
		}

		if(v.size() <= 1) continue;

		printf("infer weights with equivalent class: ");
		printv(v);
		printf("\n");

		bool b = infer_equivalent_class(v);
		if(b == true) flag = true;
	}
	return flag;
}

bool scallop2::infer_equivalent_class(const vector<int> &v)
{
	vector<int> v1;		// trustworthy vertices
	vector<int> v2;		// lest trustworthy vertices
	int sum = 0;
	bool root = false;
	int maxi = -1;
	double maxr = 0;
	for(int i = 0; i < v.size(); i++)
	{
		if(v[i] == 0 || v[i] == gr.num_vertices() - 1)
		{
			root = true;
			continue;
		}
		vertex_info vi = gr.get_vertex_info(v[i]);
		sum += vi.length;

		if(vi.reliability >= 1.0) v1.push_back(v[i]);
		if(vi.reliability >= infer_min_reliability) v2.push_back(v[i]);
		if(vi.reliability > maxr)
		{
			maxr = vi.reliability;
			maxi = v[i];
		}
	}

	//printf("INFER %lu %lu %lu max = %.2lf root = %c\n", v.size(), v1.size(), v2.size(), maxr, root ? 'T' : 'F');

	if(v1.size() == 0) v1 = v2;

	if(v1.size() == 0 && root == true && maxr >= infer_root_reliability)
	{
		v1.push_back(maxi);
	}

	if(v1.size() == 0) return false;

	double s1 = 0;
	double s2 = 0;
	for(int i = 0; i < v1.size(); i++)
	{
		double w = gr.get_vertex_weight(v1[i]);
		double d = gr.get_vertex_info(v1[i]).stddev;
		int l = gr.get_vertex_info(v1[i]).length;
		s1 += w * l / d / d;
		s2 += l / d / d;
	}
	double ave = s1 / s2;

	for(int i = 0; i < v.size(); i++)
	{
		if(v[i] == 0) continue;
		if(v[i] == gr.num_vertices() - 1) continue;
		gr.set_vertex_weight(v[i], ave);
		vertex_info vi = gr.get_vertex_info(v[i]);
		vi.infer = true;
		gr.set_vertex_info(v[i], vi);
	}
	return true;
}

bool scallop2::infer_vertices()
{
	edge_iterator it1, it2;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		bool b1 = infer_vertex_with_in_edges(i);
		bool b2 = infer_vertex_with_out_edges(i);
		if(b1 || b2) printf("infer vertex %d\n", i);
		if(b1 || b2) return true;
	}
	return false;
}

bool scallop2::infer_vertex_with_in_edges(int x)
{
	vertex_info vi = gr.get_vertex_info(x);
	if(vi.infer == true) return false;

	edge_iterator it1, it2;
	double w = 0;
	for(tie(it1, it2) = gr.in_edges(x); it1 != it2; it1++)
	{
		edge_info ei = gr.get_edge_info(*it1);
		if(ei.infer == false) return false;
		w += gr.get_edge_weight(*it1);
	}

	vi.infer = true;
	gr.set_vertex_weight(x, w);
	gr.set_vertex_info(x, vi);

	return true;
}

bool scallop2::infer_vertex_with_out_edges(int x)
{
	vertex_info vi = gr.get_vertex_info(x);
	if(vi.infer == true) return false;

	edge_iterator it1, it2;
	double w = 0;
	for(tie(it1, it2) = gr.out_edges(x); it1 != it2; it1++)
	{
		edge_info ei = gr.get_edge_info(*it1);
		if(ei.infer == false) return false;
		w += gr.get_edge_weight(*it1);
	}

	vi.infer = true;
	gr.set_vertex_weight(x, w);
	gr.set_vertex_info(x, vi);

	return true;
}

bool scallop2::infer_trivial_edges()
{
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		bool b1 = infer_trivial_in_edge(i);
		bool b2 = infer_trivial_out_edge(i);
		if(b1 || b2) printf("infer trivial edges with vertex %d\n", i);
		if(b1 || b2) return true;
	}
	return false;
}

bool scallop2::infer_trivial_in_edge(int v)
{
	vertex_info vi = gr.get_vertex_info(v);
	if(vi.infer == false) return false;
	if(gr.in_degree(v) != 1) return false;

	edge_iterator it1, it2;
	tie(it1, it2) = gr.in_edges(v);
	edge_descriptor e = (*it1);
	edge_info ei = gr.get_edge_info(e);

	if(ei.infer == true) return false;

	double w = gr.get_vertex_weight(v);

	gr.set_edge_weight(e, w);
	ei.infer = true;
	gr.set_edge_info(e, ei);

	return true;
}

bool scallop2::infer_trivial_out_edge(int v)
{
	vertex_info vi = gr.get_vertex_info(v);
	if(vi.infer == false) return false;
	if(gr.out_degree(v) != 1) return false;

	edge_iterator it1, it2;
	tie(it1, it2) = gr.out_edges(v);
	edge_descriptor e = (*it1);
	edge_info ei = gr.get_edge_info(e);

	if(ei.infer == true) return false;

	double w = gr.get_vertex_weight(v);

	gr.set_edge_weight(e, w);
	ei.infer = true;
	gr.set_edge_info(e, ei);

	return true;
}

bool scallop2::infer_edges()
{
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		bool b1 = infer_in_edges(i);
		bool b2 = infer_out_edges(i);
		if(b1 || b2) printf("infer edges with vertex %d\n", i);
		if(b1 || b2) return true;
	}
	return false;
}

bool scallop2::infer_in_edges(int v)
{
	vertex_info vi = gr.get_vertex_info(v);
	if(vi.infer == false) return false;

	double w = gr.get_vertex_weight(v);

	edge_iterator it1, it2;

	vector<edge_descriptor> vv1;	// to source s
	vector<edge_descriptor> vv2;	// to other vertices
	double ww = w;
	for(tie(it1, it2) = gr.in_edges(v); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		edge_info ei = gr.get_edge_info(e);

		if(ei.infer == true) ww -= gr.get_edge_weight(e);
		else if(e->source() == 0 && ei.length == 0) vv1.push_back(e);
		else vv2.push_back(e);
	}

	if(ww <= 0.0) return false;
	if(vv1.size() >= 1) return false;
	if(vv2.size() == 0) return false;

	distribute_weights(ww, vv2);

	return true;
}

bool scallop2::infer_out_edges(int v)
{
	vertex_info vi = gr.get_vertex_info(v);
	if(vi.infer == false) return false;

	double w = gr.get_vertex_weight(v);

	edge_iterator it1, it2;

	vector<edge_descriptor> vv1;	// to target t
	vector<edge_descriptor> vv2;	// to other vertices
	double ww = w;

	for(tie(it1, it2) = gr.out_edges(v); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		edge_info ei = gr.get_edge_info(e);

		if(ei.infer == true) ww -= gr.get_edge_weight(e);
		else if(e->target() == gr.num_vertices() - 1) vv1.push_back(e);
		else vv2.push_back(e);
	}

	if(ww <= 0.0) return false;
	if(vv1.size() >= 1) return false;
	if(vv2.size() == 0) return false;

	distribute_weights(ww, vv2);

	return true;
}

int scallop2::distribute_weights(double ww, const vector<int> &ve)
{
	VE v;
	for(int i = 0; i < ve.size(); i++)
	{
		edge_descriptor e = i2e[ve[i]];
		assert(e != null_edge);
		v.push_back(e);
	}
	distribute_weights(ww, v);
	return 0;
}

int scallop2::distribute_weights(double ww, const VE &ve)
{
	if(ve.size() == 0) return 0;

	double sum = 0;
	for(int i = 0; i < ve.size(); i++)
	{
		edge_descriptor e = ve[i];
		sum += gr.get_edge_weight(e);
	}

	if(sum <= SMIN) return 0;

	for(int i = 0; i < ve.size(); i++)
	{
		edge_descriptor e = ve[i];
		edge_info ei = gr.get_edge_info(e);
		double w2 = gr.get_edge_weight(e) / sum * ww;
		gr.set_edge_weight(e, w2);
		ei.infer = true;
		gr.set_edge_info(e, ei);
	}

	return 0;
}


int scallop2::rescale_weights()
{
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		rescale_5end_weights(i);
		rescale_3end_weights(i);
	}
	return 0;
}

int scallop2::rescale_5end_weights(int i)
{
	if(i == 0) return 0;
	if(i == gr.num_vertices() - 1) return 0;
	if(gr.get_vertex_info(i).length == 0) return 0;

	VE ss;
	gr.compute_maximum_st_path_w(ss, 0, i);
	
	int a = 0;
	for(int i = 0; i < ss.size(); i++)
	{
		edge_descriptor e = ss[i];
		int s = e->source();
		a += gr.get_edge_info(e).length;
		a += gr.get_vertex_info(s).length;
	}

	int b = average_slope_length;
	if(a >= b) return 0;

	// rescale vertex weight
	int c = a + gr.get_vertex_info(i).length;
	double w = gr.get_vertex_weight(i);
	double ww = 0;
	if(c >= b) ww = 2.0 * b * (c - a) * w / (2.0 * b * c - a * a - b * b);
	else ww = 2.0 * b * w / (a + c);

	if(gr.get_vertex_info(i).infer == false)
	{
		gr.set_vertex_weight(i, ww);
	}

	if(c >= average_slope_length) return 0;

	// rescale edge weight
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.out_edges(i); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		w = gr.get_edge_weight(e);
		ww = b * w / c;
		gr.set_edge_weight(e, ww);
	}

	return 0;
}

int scallop2::rescale_3end_weights(int i)
{
	if(i == 0) return 0;
	if(i == gr.num_vertices() - 1) return 0;
	if(gr.get_vertex_info(i).length == 0) return 0;

	VE ss;
	gr.compute_maximum_st_path_w(ss, i, gr.num_vertices() - 1);
	
	int a = 0;
	for(int i = 0; i < ss.size(); i++)
	{
		edge_descriptor e = ss[i];
		int t = e->target();
		a += gr.get_edge_info(e).length;
		a += gr.get_vertex_info(t).length;
	}

	int b = average_slope_length;
	if(a >= b) return 0;

	// rescale vertex weight
	int c = a + gr.get_vertex_info(i).length;
	double w = gr.get_vertex_weight(i);
	double ww = 0;
	if(c >= b) ww = 2.0 * b * (c - a) * w / (2.0 * b * c - a * a - b * b);
	else ww = 2.0 * b * w / (a + c);

	if(gr.get_vertex_info(i).infer == false)
	{
		gr.set_vertex_weight(i, ww);
	}

	if(c >= average_slope_length) return 0;

	// rescale edge weight
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.in_edges(i); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		w = gr.get_edge_weight(e);
		ww = b * w / c;
		gr.set_edge_weight(e, ww);
	}
	return 0;
}

bool scallop2::decompose_trivial_vertices()
{
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.in_degree(i) == 0) continue;
		bool b = decompose_trivial_vertex(i);
		if(b == true) return true;
	}
	return false;
}

bool scallop2::decompose_trivial_vertex(int i)
{
	vector<int> ve;
	return decompose_trivial_vertex(i, ve);
}

bool scallop2::decompose_trivial_vertex(int i, vector<int> &ve)
{
	if(i == 0) return false;
	if(i == gr.num_vertices() - 1) return false;
	if(gr.in_degree(i) >= 2 && gr.out_degree(i) >= 2) return false;
	if(gr.degree(i) == 0) return true;

	if(ve.size() == 0)
	{
		printf("decompose trivial vertex %d\n", i);
	}

	edge_iterator it1, it2, it3, it4;
	tie(it1, it2) = gr.in_edges(i);
	tie(it3, it4) = gr.out_edges(i);
	int e1 = e2i[*it1];
	int e2 = e2i[*it3];

	int ee = merge_adjacent_edges(e1, e2);
	ve.push_back(ee);

	return decompose_trivial_vertex(i, ve);
}

bool scallop2::remove_false_boundary_edges()
{
	edge_iterator it1, it2;
	VE v;
	for(tie(it1, it2) = gr.out_edges(0); it1 != it2; it1++)
	{
		v.push_back(*it1);
	}
	for(tie(it1, it2) = gr.in_edges(gr.num_vertices() - 1); it1 != it2; it1++)
	{
		v.push_back(*it1);
	}

	for(int i = 0; i < v.size(); i++)
	{
		edge_descriptor e = v[i];


		bool b = verify_false_boundary_edge(e);
		if(b == false) continue;

		int k = e2i[e];
		gr.remove_edge(e);
		e2i.erase(e);
		i2e[k] = null_edge;

		return true;
	}
	return false;
}

bool scallop2::verify_false_boundary_edge(edge_descriptor e)
{
	int n = gr.num_vertices() - 1;
	int k = e2i[e];
	int s = e->source();
	int t = e->target();
	int l = gr.get_edge_info(e).length;
	//printf("verify edge %d = (%d, %d), length = %d, we = %.2lf, wt = %.2lf\n", k, s, t, l, we, wt);

	if(s != 0 && t != n) return false;
	if(l != 0) return false;
	if(s == 0 && gr.in_degree(t) == 1) return false;
	if(t == n && gr.out_degree(s) == 1) return false;

	double we = gr.get_edge_weight(e);
	double ww = 0;
	if(s == 0) ww = gr.get_vertex_weight(t);
	if(t == n) ww = gr.get_vertex_weight(s);

	if(ww <= SMIN) return false;
	if(we / ww >= min_boundary_edge_weight_ratio) return false;

	printf("remove false boundary edge %d (%.2lf / %.2lf)\n", k, we, ww);
	return true;
}

int scallop2::greedy_decompose(int num)
{
	int cnt = 0;
	while(true)
	{
		if(num != -1 && cnt >= num) break;

		VE v;
		vector<int> vv;
		double w = gr.compute_maximum_path_w(v);
		if(w <= 0.0) break;
		int e = split_merge_path(v, w, vv);
		collect_path(e);
		cnt++;
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

double scallop2::compute_equation_error_ratio(equation &eqn)
{
	double w1 = 0;
	double w2 = 0;
	for(int i = 0; i < eqn.s.size(); i++)
	{
		assert(i2e[eqn.s[i]] != null_edge);
		edge_descriptor e = i2e[eqn.s[i]];
		w1 += gr.get_edge_weight(e);
	}
	for(int i = 0; i < eqn.t.size(); i++)
	{
		assert(i2e[eqn.t[i]] != null_edge);
		edge_descriptor e = i2e[eqn.t[i]];
		w2 += gr.get_edge_weight(e);
	}

	eqn.e = fabs(w1 - w2);

	if(w1 <= SMIN && w2 <= SMIN) return -1;
	double w = w1 > w2 ? w1 : w2;
	return eqn.e * 1.0 / w;
}

int scallop2::print()
{
	int n = 0;
	for(int i = 0; i < gr.num_vertices(); i++) 
	{
		if(gr.degree(i) >= 1) n++;
	}

	/*
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		printf("edge %d (%d, %d), weight = %.2lf\n", e2i[e], e->source(), e->target(), gr.get_edge_weight(e));
	}
	*/

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
		vertex_info vi = gr.get_vertex_info(i);
		double w = gr.get_vertex_weight(i);
		int l = vi.length;
		double d = vi.reliability;
		char b = vi.infer ? 'T' : 'F';
		//string s = gr.get_vertex_string(i);
		//sprintf(buf, "%d:%.0lf:%s", i, w, s.c_str());
		sprintf(buf, "%d:%.1lf:%d:%.2lf:%c", i, w, l, d, b);
		mis.insert(PIS(i, buf));
	}

	MES mes;
	for(int i = 0; i < i2e.size(); i++)
	{
		if(i2e[i] == null_edge) continue;
		double w = gr.get_edge_weight(i2e[i]);
		edge_info ei = gr.get_edge_info(i2e[i]);
		int l = ei.length;
		char b = ei.infer ? 'T' : 'F';
		sprintf(buf, "%d:%.1lf:%d:%c", i, w, l, b);
		mes.insert(PES(i2e[i], buf));
	}
	gr.draw(file, mis, mes, 4.5);
	return 0;
}
