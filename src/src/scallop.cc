#include "scallop.h"
#include "subsetsum.h"
#include "nested_graph.h"
#include "config.h"
#include "smoother.h"

#include <cstdio>
#include <iostream>
#include <cfloat>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/adjacency_list.hpp>

scallop::scallop(const string &s, splice_graph &g)
	: name(s), gr(g)
{}

scallop::~scallop()
{}

int scallop::assemble()
{
	if(algo == "scallop0") return assemble0();
	if(algo == "scallop1") return assemble1();
	if(algo == "scallop2") return assemble2();
	if(algo == "greedy") return greedy();
	return 0;
}

int scallop::assemble0()
{
	printf("\nprocess %s\n", name.c_str());
	round = 0;

	gr.draw(name + "." + tostring(round++) + ".tex");

	gr.smooth_weights();
	gr.round_weights();
	gr.remove_empty_edges();

	gr.draw(name + "." + tostring(round++) + ".tex");

	init_super_edges();
	reconstruct_splice_graph();
	gr.get_edge_indices(i2e, e2i);
	init_disjoint_sets();
	nt.build(gr);

	collect_existing_st_paths();
	print();

	printf("%s scallop0 solution %lu paths\n", name.c_str(), paths.size());

	return 0;
}

int scallop::greedy()
{
	assemble0();

	greedy_decompose();
	assert(gr.num_edges() == 0);

	print();

	printf("%s greedy solution %lu paths\n", name.c_str(), paths.size());

	return 0;
}

int scallop::assemble1()
{
	assemble0();

	iterate4();

	collect_existing_st_paths();
	print();

	printf("%s scallop1 solution %lu paths\n", name.c_str(), paths.size());

	return 0;
}

int scallop::assemble2()
{
	assemble1();

	greedy_decompose();
	assert(gr.num_edges() == 0);

	print();

	printf("%s scallop2 solution %lu paths\n", name.c_str(), paths.size());

	return 0;
}

bool scallop::iterate4()
{
	bool flag = false;
	while(true)
	{
		bool b2 = iterate3();

		vector<int> subs;
		vector<int> subt;
		bool b0 = identify_equation2(subs, subt);
		bool b1 = false;
		if(b0 == true) b1 = split_equation(subs, subt);
		if(b1 == true) 
		{
			nt.build(gr);

			printf("equation subs = (%d", subs[0]);
			for(int i = 1; i < subs.size(); i++) printf(", %d", subs[i]);
			printf("), subt = (%d", subt[0]);
			for(int i = 1; i < subt.size(); i++) printf(", %d", subt[i]);
			printf(")\n");

			if(subs.size() >= 2 || subt.size() >= 2) print();
		}

		if(b1 == true || b2 == true) flag = true;
		if(b1 == false && b2 == false) break;
	}
	return flag;
}

bool scallop::iterate3()
{
	bool flag = false;
	while(true)
	{
		bool b2 = iterate2();

		int ex, ey;
		bool b1 = compute_closest_equal_edges(ex, ey);

		if(b1 == true)
		{
			printf("shortest equal path (%d, %d)\n", ex, ey);
			connect_equal_edges(ex, ey);
			nt.build(gr);
			print();
		}

		if(b1 == true || b2 == true) flag = true;
		if(b1 == false && b2 == false) break;
	}
	return flag;
}

bool scallop::iterate2()
{
	bool flag = false;
	while(true)
	{
		bool b2 = iterate1();

		vector<int> subs;
		vector<int> subt;
		bool b1 = identify_equation1(subs, subt);
		if(b1 == true) 
		{
			split_equation(subs, subt);
			nt.build(gr);

			if(subs.size() >= 2 || subt.size() >= 2)
			{
				printf("equation subs = (%d", subs[0]);
				for(int i = 1; i < subs.size(); i++) printf(", %d", subs[i]);
				printf("), subt = (%d", subt[0]);
				for(int i = 1; i < subt.size(); i++) printf(", %d", subt[i]);
				printf(")\n");

				print();
			}
		}

		if(b1 == true || b2 == true) flag = true;
		if(b1 == false && b2 == false) break;
	}
	return flag;
}

bool scallop::iterate1()
{
	bool flag = false;
	while(true)
	{
		bool b2 = decompose_trivial_vertices();
		nt.build(gr);
		if(b2 == true)
		{
			print();
		}

		int ex, ey;
		vector<PI> p;
		bool b1 = identify_linkable_edges(ex, ey, p);
		if(b1 == true)
		{
			printf("link edges = (%d, %d) using path = (#", ex, ey);
			for(int i = 0; i < p.size(); i++) printf(", (%d,%d)", p[i].first, p[i].second);
			printf(")\n");

			build_adjacent_equal_edges(p);
			connect_adjacent_equal_edges(ex, ey);
			nt.build(gr);

			print();
		}

		if(b1 == true || b2 == true) flag = true;
		if(b1 == false && b2 == false) break;
	}
	return flag;
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

int scallop::init_disjoint_sets()
{
	ds = disjoint_sets_t(gr.num_edges() * gr.num_vertices());
	for(int i = 0; i < gr.num_edges(); i++)
	{
		ds.make_set(i);
	}
	return 0;
}

vector<int> scallop::compute_representatives()
{
	vector<int> v;
	vector< vector<int> > vv = get_disjoint_sets(ds, i2e.size());
	for(int i = 0; i < vv.size(); i++)
	{
		if(vv[i].size() <= 0) continue;
		int k = -1;
		for(int j = 0; j < vv[i].size(); j++)
		{
			int e = vv[i][j];
			if(i2e[e] == null_edge) continue;
			k = e;
			break;
		}
		if(k == -1) continue;
		v.push_back(k);
	}
	return v;
}

vector< vector<int> > scallop::compute_disjoint_sets()
{
	vector< vector<int> > xx;
	vector< vector<int> > vv = get_disjoint_sets(ds, i2e.size());
	for(int i = 0; i < vv.size(); i++)
	{
		if(vv[i].size() == 0) continue;
		vector<int> v;
		for(int j = 0; j < vv[i].size(); j++)
		{
			int e = vv[i][j];
			if(i2e[e] == null_edge) continue;
			v.push_back(e);
		}
		if(v.size() <= 0) continue;
		xx.push_back(v);
	}
	return xx;
}

set<int> scallop::compute_singletons()
{
	set<int> s;
	vector< vector<int> > vv = compute_disjoint_sets();
	for(int i = 0; i < vv.size(); i++)
	{
		assert(vv[i].size() >= 1);
		if(vv[i].size() >= 2) continue;
		s.insert(vv[i][0]);
	}
	return s;
}

int scallop::connect_equal_edges(int x, int y)
{
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

	VE p;
	int l = gr.compute_shortest_path_w(xt, ys, wx, p);
	assert(l >= 0);

	p.insert(p.begin(), xx);
	p.insert(p.end(), yy);

	return connect_path(p, wx);
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

	ds.make_set(n);
	ds.union_set(n, x);
	ds.union_set(n, y);

	e2i.erase(xx);
	e2i.erase(yy);
	i2e[x] = null_edge;
	i2e[y] = null_edge;
	gr.remove_edge(xx);
	gr.remove_edge(yy);

	return n;
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

	int ee = split_edge(exi, wy);

	if(ee == exi)
	{
		ds.union_set(exi, eyi);
		return exi;
	}

	ds.make_set(ee);
	ds.make_set(ee + 1);
	ds.union_set(ee, eyi);

	return ee;
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

bool scallop::split_equation(const vector<int> &subs, const vector<int> &subt)
{
	bool b = split_equation_maxflow(subs, subt);
	if(b == false) return false;
	split_equation_greedy(subs, subt);
	return true;
}

int scallop::split_equation_greedy(const vector<int> &subs, const vector<int> &subt)
{
	vector<int> vv;
	for(int i = 0; i < subs.size(); i++) vv.push_back(subs[i]);
	for(int i = 0; i < subt.size(); i++) vv.push_back(subt[i]);

	directed_graph g;
	for(int i = 0; i < vv.size(); i++) g.add_vertex();

	for(int i = 0; i < subs.size(); i++)
	{
		for(int j = subs.size(); j < vv.size(); j++)
		{
			bool f1 = gr.check_path(i2e[vv[i]], i2e[vv[j]]);
			bool f2 = gr.check_path(i2e[vv[j]], i2e[vv[i]]);
			assert(f1 == false || f2 == false);
			if(f1 == false && f2 == false) continue;
			g.add_edge(i, j);
		}
	}

	while(g.num_edges() >= 1)
	{
		bool flag = false;
		for(int i = 0; i < subs.size(); i++)
		{
			if(g.out_degree(i) == 0) continue;
			if(g.out_degree(i) >= 2) continue;
			edge_iterator it1, it2;
			tie(it1, it2) = g.out_edges(i);
			int j = (*it1)->target();
			double wi = gr.get_edge_weight(i2e[vv[i]]);
			double wj = gr.get_edge_weight(i2e[vv[j]]);
			//printf(" A: edge %d:%d:%.0lf -> %d:%d:%.0lf\n", i, vv[i], wi, j, vv[j], wj);
			assert(wi <= wj);
			int e = split_edge(vv[j], vv[i]);
			g.clear_vertex(i);
			if(e == vv[j]) assert(g.degree(j) == 0);
			vv[j] = e + 1;
			flag = true;
			break;
		}
		if(flag == true) continue;

		for(int j = subs.size(); j < vv.size(); j++)
		{
			if(g.in_degree(j) == 0) continue;
			if(g.in_degree(j) >= 2) continue;
			edge_iterator it1, it2;
			tie(it1, it2) = g.in_edges(j);
			int i = (*it1)->source();
			double wi = gr.get_edge_weight(i2e[vv[i]]);
			double wj = gr.get_edge_weight(i2e[vv[j]]);
			//printf(" B: edge %d:%d:%.0lf -> %d:%d:%.0lf\n", i, vv[i], wi, j, vv[j], wj);
			assert(wi >= wj);
			int e = split_edge(vv[i], vv[j]);
			g.clear_vertex(j);
			if(e == vv[i]) assert(g.degree(i) == 0);
			vv[i] = e + 1;
			flag = true;
			break;
		}
		if(flag == true) continue;

		edge_iterator it1, it2;
		tie(it1, it2) = g.edges();
		int i = (*it1)->source();
		int j = (*it1)->target();
		assert(i < subs.size());
		assert(j >= subs.size());
		double wi = gr.get_edge_weight(i2e[vv[i]]);
		double wj = gr.get_edge_weight(i2e[vv[j]]);
		//printf(" C: edge %d:%d:%.0lf -> %d:%d:%.0lf\n", i, vv[i], wi, j, vv[j], wj);
		if(wi >= wj)
		{
			int e = split_edge(vv[i], vv[j]);
			g.clear_vertex(j);
			if(e == vv[i]) assert(g.degree(i) == 0);
			vv[i] = e + 1;
		}
		else
		{
			int e = split_edge(vv[j], vv[i]);
			g.clear_vertex(i);
			if(e == vv[j]) assert(g.degree(j) == 0);
			vv[j] = e + 1;
		}
	}
	return 0;
}

bool scallop::split_equation_maxflow(const vector<int> &subs, const vector<int> &subt)
{
	using namespace boost;

	typedef adjacency_list_traits<vecS, vecS, bidirectionalS> Traits;
	typedef adjacency_list<vecS, vecS, bidirectionalS, 
			property<vertex_name_t, std::string>,
			property<edge_capacity_t, int,
			property<edge_residual_capacity_t, int,
			property<edge_reverse_t, Traits::edge_descriptor> > > > Graph;

	Graph g;

	property_map<Graph, edge_capacity_t>::type capacity = get(edge_capacity, g);
	property_map<Graph, edge_reverse_t>::type rev = get(edge_reverse, g);
	property_map<Graph, edge_residual_capacity_t>::type residual_capacity = get(edge_residual_capacity, g);

	// vertices
	add_vertex(g);
	for(int i = 0; i < subs.size(); i++) add_vertex(g);
	for(int i = 0; i < subt.size(); i++) add_vertex(g);
	add_vertex(g);

	int n = num_vertices(g);

	int s = 0;
	int t = n - 1;

	// edges from s to subs
	Traits::edge_descriptor e1, e2;
	bool b1, b2;
	int sums = 0;
	int sumt = 0;
	for(int i = 0; i < subs.size(); i++)
	{
		tie(e1, b1) = add_edge(s, i + 1, g);
		tie(e2, b2) = add_edge(i + 1, s, g);
		int w = (int)(gr.get_edge_weight(i2e[subs[i]]));
		capacity[e1] = w;
		capacity[e2] = 0;
		rev[e1] = e2;
		rev[e2] = e1;
		sums += w;
	}

	// edges from subt to t
	for(int i = 0; i < subt.size(); i++)
	{
		tie(e1, b1) = add_edge(i + 1 + subs.size(), t, g);
		tie(e2, b2) = add_edge(t, i + 1 + subs.size(), g);
		int w = (int)(gr.get_edge_weight(i2e[subt[i]]));
		capacity[e1] = w;
		capacity[e2] = 0;
		rev[e1] = e2;
		rev[e2] = e1;
		sumt += w;
	}

	// edges between subs and subt
	for(int i = 0; i < subs.size(); i++)
	{
		for(int j = 0; j < subt.size(); j++)
		{
			bool f1 = gr.check_path(i2e[subs[i]], i2e[subt[j]]);
			bool f2 = gr.check_path(i2e[subt[j]], i2e[subs[i]]);
			assert(f1 == false || f2 == false);
			if(f1 == false && f2 == false) continue;

			tie(e1, b1) = add_edge(i + 1, j + 1 + subs.size(), g);
			tie(e2, b2) = add_edge(j + 1 + subs.size(), i + 1, g);
			capacity[e1] = sums + sumt;
			capacity[e2] = 0;
			rev[e1] = e2;
			rev[e2] = e1;
		}
	}

	int sum = (sums < sumt) ? sums : sumt;

	int flow = push_relabel_max_flow(g, s, t);

	assert(flow <= sum);
	if(flow != sum) return false;
	else return true;


	// print flow
	Traits::edge_descriptor ee;
	bool bb;
	for(int i = 0; i < subs.size(); i++)
	{
		tie(ee, bb) = edge(s, i + 1, g);
		printf("edge %3d -> %3d, capacity = %5d, flow = %5d\n", 
				0, i + 1, capacity[ee], capacity[ee] - residual_capacity[ee]);
	}
	for(int i = 0; i < subt.size(); i++)
	{
		tie(ee, bb) = edge(i + 1 + subs.size(), t, g);
		printf("edge %3lu -> %3d, capacity = %5d, flow = %5d\n", 
				i + 1 + subs.size(), t, capacity[ee], capacity[ee] - residual_capacity[ee]);
	}
	for(int i = 0; i < subs.size(); i++)
	{
		for(int j = 0; j < subt.size(); j++)
		{
			tie(ee, bb) = edge(i + 1, j + 1 + subs.size(), g);
			if(bb == false) continue;
			printf("edge %3d -> %3lu, capacity = %5d, flow = %5d\n", 
					i + 1, j + 1 + subs.size(), capacity[ee], capacity[ee] - residual_capacity[ee]);
		}
	}


	vector< vector<int> > st(subs.size());
	vector< vector<int> > ts(subt.size());
	graph_traits<Graph>::out_edge_iterator oi1, oi2;
	for(int i = 0; i < subs.size(); i++) 
	{
		st[i].assign(subt.size(), -1);

		int ee = subs[i];
		for (tie(oi1, oi2) = out_edges(i + 1, g); oi1 != oi2; oi1++)
		{
			int ss = source(*oi1, g) - 1;
			int tt = target(*oi1, g) - 1 - subs.size();
			assert(ss == i);
			if(tt < 0) continue;
			int f = capacity[*oi1] - residual_capacity[*oi1];
			assert(f >= 0);
			if(f == 0) continue;
			st[ss][tt] = split_edge(ee, (double)(f));
			ee = st[ss][tt] + 1;
		}
	}
	graph_traits<Graph>::in_edge_iterator ii1, ii2;
	for(int i = 0; i < subt.size(); i++)
	{
		ts[i].assign(subs.size(), -1);

		int ee = subt[i];
		for (tie(ii1, ii2) = in_edges(i + 1 + subs.size(), g); ii1 != ii2; ii1++)
		{
			int ss = source(*ii1, g) - 1;
			int tt = target(*ii1, g) - 1 - subs.size();
			assert(tt == i);
			if(ss == t - 1) continue;
			int f = capacity[*ii1] - residual_capacity[*ii1];
			assert(f >= 0);
			if(f == 0) continue;
			ts[tt][ss] = split_edge(ee, (double)(f));
			ee = ts[tt][ss] + 1;
		}
	}

	for(int i = 0; i < subs.size(); i++)
	{
		for(int j = 0; j < subt.size(); j++)
		{
			if(st[i][j] == -1) assert(ts[j][i] == -1);
			if(ts[j][i] == -1) assert(st[i][j] == -1);
			if(st[i][j] == -1) continue;
			printf("connect %d and %d\n", st[i][j], ts[j][i]);
			split_edge(st[i][j], ts[j][i]);
		}
	}
	return true;
}

bool scallop::identify_equation1(vector<int> &subs, vector<int> &subt)
{
	vector<PI> p;
	for(int i = 0; i < i2e.size(); i++)
	{
		if(i2e[i] == null_edge) continue;
		int w = (int)(gr.get_edge_weight(i2e[i]));
		p.push_back(PI(w, i));
	}
	sort(p.begin(), p.end());

	subs.clear();
	subt.clear();
	for(int i = 0; i < p.size(); i++)
	{
		int e = p[i].second;
		assert(i2e[e] != null_edge);
		vector<int> s;
		s.push_back(e);
		vector<int> t;
		bool b = identify_equation(s, t);
		if(b == false) continue;
		if(subt.size() == 0 || t.size() < subt.size())
		{
			subs = s;
			subt = t;
		}
	}

	if(subs.size() >= 1) return true;
	else return false;
}

bool scallop::identify_equation2(vector<int> &subs, vector<int> &subt)
{
	vector<PI> p;
	for(int i = 0; i < i2e.size(); i++)
	{
		if(i2e[i] == null_edge) continue;
		int w = (int)(gr.get_edge_weight(i2e[i]));
		p.push_back(PI(w, i));
	}
	sort(p.begin(), p.end());

	subs.clear();
	subt.clear();
	for(int i = 0; i < p.size(); i++)
	{
		for(int j = i + 1; j < p.size(); j++)
		{
			int e1 = p[i].second;
			int e2 = p[j].second;
			if(i2e[e1]->source() == i2e[e2]->source() && gr.out_degree(i2e[e1]->source()) == 2) continue;
			if(i2e[e1]->target() == i2e[e2]->target() && gr.in_degree(i2e[e1]->target()) == 2) continue;
			vector<int> s;
			s.push_back(e1);
			s.push_back(e2);
			vector<int> t;
			bool b = identify_equation(s, t);
			if(b == false) continue;
			if(subt.size() == 0 || t.size() < subt.size())
			{
				subs = s;
				subt = t;
			}
		}
	}

	if(subs.size() >= 1) return true;
	else return false;
}

bool scallop::identify_equation(const vector<int> &subs, vector<int> &subt)
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

	set<int> sr = compute_singletons();

	set<int> ss(subs.begin(), subs.end());
	vector<PI> xi;
	for(SE::iterator it = ff.begin(); it != ff.end(); it++)
	{
		if(sr.find(e2i[*it]) == sr.end()) continue;
		if(ss.find(e2i[*it]) != ss.end()) continue;
		int ww = (int)(gr.get_edge_weight(*it));
		if(ww > sw) continue;
		xi.push_back(PI(ww, e2i[*it]));
	}

	sort(xi.begin(), xi.end());

	xi.push_back(PI(sw, e));

	vector<int> v;
	for(int i = 0; i < xi.size(); i++)
	{
		v.push_back(xi[i].first);
	}

	subsetsum sss(v);
	sss.solve();

	subt.clear();
	for(int j = 0; j < sss.subset.size(); j++)
	{
		int k = sss.subset[j];
		assert(xi[k].second != -1);
		subt.push_back(xi[k].second);
	}

	int err = (int)fabs(sss.opt - sw);
	if(err >= 1) return false;

	assert(subt.size() >= 1);

	/*
	printf("closest subset error = %d, S = (", err);
	for(int i = 0; i < subs.size() - 1; i++) printf("%d:%.0lf, ", subs[i], gr.get_edge_weight(i2e[subs[i]]));
	printf("%d:%.0lf), T = (", subs[subs.size() - 1], gr.get_edge_weight(i2e[subs[subs.size() - 1]]));
	for(int i = 0; i < subt.size() - 1; i++) printf("%d:%.0lf, ", subt[i], gr.get_edge_weight(i2e[subt[i]]));
	printf("%d:%.0lf)\n", subt[subt.size() - 1], gr.get_edge_weight(i2e[subt[subt.size() - 1]]));
	*/

	return true;
}

bool scallop::check_linkable(int ex, int ey, vector<PI> &p)
{
	assert(i2e[ex] != null_edge);
	assert(i2e[ey] != null_edge);

	int xs = i2e[ex]->source();
	int xt = i2e[ex]->target();
	int ys = i2e[ey]->source();
	int yt = i2e[ey]->target();

	vector<PI> xp, yp;
	bool b = nt.link(xs, xt, ys, yt, xp, yp);
	if(b == false) return false;

	p = xp;
	p.insert(p.end(), yp.begin(), yp.end());

	return true;
}

bool scallop::identify_linkable_edges(int &ex, int &ey, vector<PI> &p)
{
	ex = ey = -1;
	p.clear();
	vector< vector<int> > vv = compute_disjoint_sets();
	bool flag = false;
	for(int i = 0; i < vv.size(); i++)
	{
		vector<int> &v = vv[i];
		if(v.size() == 1) continue;
		for(int j = 0; j < v.size(); j++)
		{
			for(int k = j + 1; k < v.size(); k++)
			{
				bool b = check_linkable(v[j], v[k], p);
				//printf(" check %d and %d = %c\n", v[j], v[k], b ? 'T' : 'F');
				if(b == false) continue;
				ex = v[j];
				ey = v[k];
				flag = true;
				break;
			}
			if(flag == true) break;
		}
		if(flag == true) break;
	}
	if(flag == false) return false;
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
			printf("decompose trivial vertex %d\n", i);

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
			printf("decompose trivial vertex %d\n", i);

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

bool scallop::compute_closest_equal_edges(int &ex, int &ey)
{
	ex = ey = -1;
	vector< vector<int> > vv = compute_disjoint_sets();
	int min = INT_MAX;
	for(int i = 0; i < vv.size(); i++)
	{
		vector<int> &v = vv[i];
		if(v.size() == 1) continue;
		for(int j = 0; j < v.size(); j++)
		{
			int xs = i2e[v[j]]->source();
			int xt = i2e[v[j]]->target();
			double wx = gr.get_edge_weight(i2e[v[j]]);
			for(int k = j + 1; k < v.size(); k++)
			{
				int ys = i2e[v[k]]->source();
				int yt = i2e[v[k]]->target();
				double wy = gr.get_edge_weight(i2e[v[k]]);

				assert(fabs(wx - wy) <= SMIN);

				int pxy = gr.compute_shortest_path_w(xt, ys, wx);
				int pyx = gr.compute_shortest_path_w(yt, xs, wx);
				
				assert(pxy == -1 || pyx == -1);
				if(pxy == -1 && pyx == -1) continue;

				if(pxy >= 0 && pxy < min)
				{
					min = pxy;
					ex = v[j];
					ey = v[k];
				}
				else if(pyx >= 0 && pyx < min)
				{
					min = pyx;
					ex = v[k];
					ey = v[j];
				}
			}
		}
	}
	if(min == INT_MAX) return false;
	else return true;
}

int scallop::greedy_decompose()
{
	while(true)
	{
		VE v;
		double w = gr.compute_maximum_path_w(v);
		if(w <= 0.5) break;
		int e = connect_path(v, w);
		collect_path(e);
	}
	return 0;
}

bool scallop::identify_optimal_paths()
{
	bool flag = false;
	while(true)
	{
		VE v;
		bool b = gr.compute_optimal_path(v);
		if(b == false) break;
		flag = true;

		double w = gr.compute_minimum_weight(v);
		int e = connect_path(v, w);

		printf("connect optimal path %d, weight = %.2lf\n", e, w);
		print();

		collect_path(e);
	}
	return flag;
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

int scallop::print()
{
	/*
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		vector<int> v = mev[*it1];
		printf("vertices in edge %d = ", e2i[*it1]);
		for(int i = 0; i < v.size(); i++) printf("%d ", v[i]);
		printf("\n");
	}

	vector< vector<int> > vv = compute_disjoint_sets();
	for(int i = 0; i < vv.size(); i++)
	{
		vector<int> v = vv[i];
		assert(v.size() >= 1);

		if(v.size() == 1) continue;
		int w = (int)(gr.get_edge_weight(i2e[v[0]]));

		printf("edge set %d, weight = %d, #edges = %lu, set = (%d", i, w, v.size(), v[0]);
		for(int j = 1; j < v.size(); j++) printf(", %d", v[j]);
		printf(")\n");
	}
	*/
	int n = 0;
	for(int i = 0; i < gr.num_vertices(); i++) 
	{
		if(gr.degree(i) >= 1) n++;
	}
	printf("statistics: %lu edges, %d vertices\n", gr.num_edges(), n);
	printf("finish round %d\n\n", round);

	if(output_tex_files == true)
	{
		char buf[1024];
		sprintf(buf, "%s.gr.%d.tex", name.c_str(), round);
		draw_splice_graph(buf);

		/*
		sprintf(buf, "%s.nt.%d.tex", name.c_str(), round);
		nt.draw(buf);
		*/
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
