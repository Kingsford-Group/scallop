#include "scallop3.h"
#include "config.h"
#include "gurobi_c++.h"
#include "binomial.h"
#include "smoother.h"

#include <cstdio>
#include <iostream>
#include <cfloat>
#include <algorithm>

scallop3::scallop3()
{}

scallop3::scallop3(const string &s, const splice_graph &g, const hyper_set &h)
	: name(s), gr(g), hs(h)
{
	round = 0;
	if(output_tex_files == true) gr.draw(name + "." + tostring(round++) + ".tex");
	gr.get_edge_indices(i2e, e2i);
	add_pseudo_hyper_edges();
	hs.build(gr, e2i);
	init_super_edges();
	init_vertex_map();
	init_inner_weights();
	print();
}

scallop3::~scallop3()
{
}

int scallop3::assemble()
{
	classify();

	while(true)
	{
		bool b	= false;

		/*
		b = resolve_ignorable_edges();
		if(b == true) print();
		if(b == true) continue;
		*/

		b = resolve_nontrivial_vertex(true, true);
		if(b == true) print();
		if(b == true) continue;

		b = resolve_hyper_tree();
		if(b == true) print();
		if(b == true) continue;

		b = resolve_trivial_vertex(true);
		if(b == true) print();
		if(b == true) continue;

		b = resolve_nontrivial_vertex(true, false);
		if(b == true) print();
		if(b == true) continue;

		/*
		b = resolve_trivial_vertex(false);
		if(b == true) print();
		if(b == true) continue;
		*/

		b = resolve_nontrivial_vertex(false, false);
		if(b == true) print();
		if(b == true) continue;

		b = resolve_hyper_edge1();
		if(b == true) print();
		if(b == true) continue;

		b = resolve_hyper_edge0();
		if(b == true) print();
		if(b == true) continue;

		break;
	}

	collect_existing_st_paths();

	print();

	greedy_decompose(-1);

	return 0;
}

int scallop3::assert_weights()
{
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		assert(w >= 0.99);
	}
	return 0;
}

int scallop3::refine_splice_graph()
{
	while(true)
	{
		bool b = false;
		edge_iterator it1, it2;
		for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
		{
			int s = (*it1)->source();
			int t = (*it1)->target();
			int e = e2i[*it1];
			if(s == 0) continue;
			if(t == gr.num_vertices() - 1) continue;
			if(gr.in_degree(s) >= 1 && gr.out_degree(t) >= 1) continue;

			printf("refine graph by removing edge %d\n", e);
			assert(false);

			remove_edge(e);
			hs.remove(e);
			b = true;
			break;
		}
		if(b == false) break;
	}
	return 0;
}

bool scallop3::resolve_nontrivial_vertex(bool split, bool hyper)
{
	int root = -1;
	double ratio1 = -1;
	vector<equation> eqns;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.in_degree(i) <= 1) continue;
		if(gr.out_degree(i) <= 1) continue;

		vector<PI> p = hs.get_routes(i, gr, e2i);

		if(split == true && hyper == true && p.size() <= 0) continue;
		if(split == true && hyper == false&& p.size() >= 1) continue;

		router rt(i, gr, e2i, i2e, p);
		rt.build();
		//rt.print();

		if(rt.status != 4) continue;
		assert(rt.ratio >= 0);
		assert(rt.eqns.size() == 2);

		if(ratio1 >= 0 && ratio1 < rt.ratio) continue;

		root = i;
		ratio1 = rt.ratio;
		eqns = rt.eqns;
	}

	if(root == -1) return false;

	int se;
	double ratio2 = compute_smallest_edge(root, se);
	double sw = gr.get_edge_weight(i2e[se]);

	double ratio = (ratio1 < ratio2) ? ratio1 : ratio2;
	if(ratio > max_split_error_ratio) return false;

	if(split == true && ratio1 <= ratio2)
	{
		printf("split %s vertex %d, ratio = %.2lf / %.2lf, degree = (%d, %d)\n", 
				hyper ? "hyper" : "normal", root, ratio1, ratio2, gr.in_degree(root), gr.out_degree(root));

		for(int i = 0; i < eqns.size(); i++) eqns[i].print(99);

		equation &eqn = eqns[0];
		assert(eqn.s.size() >= 1);
		assert(eqn.t.size() >= 1);

		split_vertex(root, eqn.s, eqn.t);

		return true;
	}

	bool b = true;
	if(hs.left_extend(se) && hs.right_extend(se)) b = false;
	if(gr.in_degree(i2e[se]->target()) <= 1) b = false;
	if(gr.out_degree(i2e[se]->source()) <= 1) b = false;
	if(split == false && ratio2 <= ratio1 && b == true)
	{
		assert(se >= 0);
		printf("remove %s small edge %d of vertex %d, weight = %.2lf, ratio = %.2lf / %.2lf, degree = (%d, %d)\n", 
				hyper ? "hyper" : "normal", se, root, sw, ratio1, ratio2, gr.in_degree(root), gr.out_degree(root));

		remove_edge(se);
		hs.remove(se);

		return true;
	}

	return false;
}

bool scallop3::resolve_hyper_tree()
{
	int root = -1;
	undirected_graph ug;
	vector<int> u2e;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.in_degree(i) <= 1) continue;
		if(gr.out_degree(i) <= 1) continue;

		vector<PI> p = hs.get_routes(i, gr, e2i);
		router rt(i, gr, e2i, i2e, p);
		rt.build();

		if(rt.status != 1) continue;
		balance_vertex(i);
		bool b = balance_vertex(rt.ug, rt.u2e);
		assert(b == true);

		root = i;
		ug = rt.ug;
		u2e = rt.u2e;
		break;
	}

	if(root == -1) return false;

	printf("resolve hyper tree %d, degree = (%d, %d)\n", root, gr.in_degree(root), gr.out_degree(root));

	decompose_tree(ug, u2e);
	assert(gr.degree(root) == 0);
	return true;
}

bool scallop3::resolve_hyper_edge0()
{
	int ee1 = -1, ee2 = -1, root = -1;
	for(int i = 1; i < gr.num_vertices(); i++)
	{
		ee1 = ee2 = root = -1;
		double ww = 0;
		edge_iterator it1, it2;
		for(tie(it1, it2) = gr.in_edges(i); it1 != it2; it1++)
		{
			int e1 = e2i[*it1];
			int e2 = -1;
			double w1 = gr.get_edge_weight(*it1);
			double w2 = 0;
			set<int> s = hs.get_successors(e1);
			if(s.size() <= 0) continue;
			for(set<int>::iterator it = s.begin(); it != s.end(); it++)
			{
				double w = gr.get_edge_weight(i2e[*it]);
				if(hs.left_extend(e1)) continue;
				if(hs.right_extend(*it)) continue;
				if(w <= w2) continue;
				w2 = w;
				e2 = (*it);
			}
			if(e1 == -1 || e2 == -1) continue;
			if(w1 <= ww || w2 <= ww) continue;
			ee1 = e1;
			ee2 = e2;
			ww = (w1 < w2) ? w1 : w2;
		}
		if(ee1 == -1 || ee2 == -1) continue;
		root = i;
		break;
	}

	if(root == -1) return false;

	balance_vertex(root);

	double ww1 = gr.get_edge_weight(i2e[ee1]);
	double ww2 = gr.get_edge_weight(i2e[ee2]);
	double ww = (ww1 <= ww2) ? ww1 : ww2;

	if(ww1 <= ww2) assert(hs.left_extend(ee1) == false);
	if(ww2 <= ww1) assert(hs.right_extend(ee2) == false);

	int k1 = split_edge(ee1, ww);
	int k2 = split_edge(ee2, ww);
	int x = merge_adjacent_equal_edges(k1, k2);

	printf(" resolve hyper edge (%d, %d) of vertex %d, weight = (%.2lf, %.2lf) -> (%d, %d) -> %d\n", ee1, ee2, root, ww1, ww2, k1, k2, x);

	hs.replace(ee1, ee2, x);
	if(k1 == ee1) hs.remove(ee1);
	if(k2 == ee2) hs.remove(ee2);

	return true;
}

bool scallop3::resolve_hyper_edge1()
{
	edge_iterator it1, it2;
	vector<int> v1, v2;
	int root = -1;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		int e = e2i[*it1];
		set<int> s;

		s = hs.get_successors(e);
		if(s.size() >= 2 && hs.right_extend(s) == false)
		{
			v1.push_back(e);
			v2.insert(v2.begin(), s.begin(), s.end());
			root = (*it1)->target();
			break;
		}

		s = hs.get_predecessors(e);
		if(s.size() >= 2 && hs.left_extend(s) == false)
		{
			v1.insert(v1.begin(), s.begin(), s.end());
			v2.push_back(e);
			root = (*it1)->source();
			break;
		}
	}

	if(v1.size() == 0 || v2.size() == 0) return false;

	printf("resolve hyper edge ( ");
	printv(v1);
	printf("), ( ");
	printv(v2);
	printf(")\n");

	assert(v1.size() == 1 || v2.size() == 1);

	balance_vertex(root);

	vector<double> w1;
	vector<double> w2;
	double sum1 = 0, sum2 = 0;
	for(int i = 0; i < v1.size(); i++)
	{
		double w = gr.get_edge_weight(i2e[v1[i]]);
		w1.push_back(w);
		sum1 += w;
	}
	for(int i = 0; i < v2.size(); i++)
	{
		double w = gr.get_edge_weight(i2e[v2[i]]);
		w2.push_back(w);
		sum2 += w;
	}

	double sum = 0.5 * (gr.get_in_weights(root) + gr.get_out_weights(root));

	double r1 = (sum1 < sum2) ? 1.0 : sum2 / sum1;
	double r2 = (sum1 > sum2) ? 1.0 : sum1 / sum2;

	for(int i = 0; i < w1.size(); i++) w1[i] *= r1;
	for(int i = 0; i < w2.size(); i++) w2[i] *= r2;

	set<int> ss;
	for(int i = 0; i < w1.size(); i++)
	{
		for(int j = 0; j < w2.size(); j++)
		{
			double w = (w1[i] < w2[j]) ? w1[i] : w2[j];

			double t1 = gr.get_edge_weight(i2e[v1[i]]);
			double t2 = gr.get_edge_weight(i2e[v2[j]]);
			int k1 = split_edge(v1[i], w);
			int k2 = split_edge(v2[j], w);
			int x = merge_adjacent_equal_edges(k1, k2);

			printf(" split (%d, %d), w = %.2lf, weight = (%.2lf, %.2lf), (%.2lf, %.2lf) -> (%d, %d) -> %d\n", v1[i], v2[j], 
					w, w1[i], w2[j], t1, t2, k1, k2, x);

			hs.replace(v1[i], v2[j], x);
			if(k1 == v1[i]) hs.remove(v1[i]);
			if(k2 == v2[j]) hs.remove(v2[j]);
			//if(k1 == v1[i]) hs.replace(v1[i], x);
			//if(k2 == v2[j]) hs.replace(v2[j], x);
		}
	}

	return true;
}

bool scallop3::resolve_ignorable_edges()
{
	bool flag = false;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.in_degree(i) <= 1) continue;
		if(gr.out_degree(i) <= 1) continue;

		int ei;
		double ratio = compute_smallest_edge(i, ei);
		edge_descriptor e = i2e[ei];
		int s = e->source();
		int t = e->target();

		//if(ratio > max_split_error_ratio) continue;

		double w = gr.get_edge_weight(e);
		if(w > max_ignorable_edge_weight) continue;
		if(hs.left_extend(ei) == true && hs.right_extend(ei) == true) continue;
		if(s == i && gr.in_degree(t) <= 1) continue;
		if(t == i && gr.out_degree(s) <= 1) continue;

		printf("remove ignorable edge %d of vertex %d, weight = %.2lf, ratio = %.2lf, degree = (%d, %d)\n", 
				ei, i, w, ratio, gr.in_degree(i), gr.out_degree(i));

		remove_edge(ei);
		hs.remove(ei);
		flag = true;
	}
	return flag;
}

bool scallop3::resolve_trivial_vertex(bool split)
{
	int root = -1;
	double ratio = -1;
	int se = -1;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.degree(i) == 0) continue;
		if(gr.in_degree(i) >= 2 && gr.out_degree(i) >= 2) continue;

		int e;
		double r = compute_smallest_edge(i, e);
		if(ratio >= 0 && ratio < r) continue;

		root = i;
		ratio = r;
		se = e;
	}

	if(root == -1) return false;

	if(split == true)
	{
		printf("resolve trivial vertex %d, ratio = %.2lf, degree = (%d, %d)\n", root, ratio, gr.in_degree(root), gr.out_degree(root));

		decompose_trivial_vertex(root);
		assert(gr.degree(root) == 0);
		return true;
	}

	bool b = true;
	double ww = gr.get_edge_weight(i2e[se]);
	if(hs.left_extend(se) && hs.right_extend(se)) b = false;
	if(gr.in_degree(i2e[se]->target()) <= 1) b = false;
	if(gr.out_degree(i2e[se]->source()) <= 1) b = false;
	if(ww > max_removable_edge_weight) b = false;
	if(ratio > max_split_error_ratio) b = false;

	if(b == true && split == false)
	{
		printf("remove trivial small edge %d of vertex %d, weight = %.2lf, ratio = %.2lf, degree = (%d, %d)\n", 
				se, root, ww, ratio, gr.in_degree(root), gr.out_degree(root));

		remove_edge(se);
		hs.remove(se);
		return true;
	}

	return false;
}

int scallop3::classify()
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

	//assert(p0 >= p1);

	bool b = (p0 <= p1) ? true : false;

	printf("\nprocess %s %s\n", name.c_str(), b ? "TRIVIAL" : "NORMAL");

	if(p0 == p1) return TRIVIAL;
	else return NORMAL;
}

int scallop3::add_pseudo_hyper_edges()
{
	for(int k = 1; k < gr.num_vertices() - 1; k++)
	{
		int s = -1, t = -1;
		double w1 = 0, w2 = 0;
		edge_iterator it1, it2;
		for(tie(it1, it2) = gr.in_edges(k); it1 != it2; it1++)
		{
			double w = gr.get_edge_weight(*it1);
			if(w <= w1) continue;
			w1 = w;
			s = (*it1)->source();
		}
		for(tie(it1, it2) = gr.out_edges(k); it1 != it2; it1++)
		{
			double w = gr.get_edge_weight(*it1);
			if(w <= w2) continue;
			w2 = w;
			t = (*it1)->target();
		}
		if(s == -1 || t == -1) continue;
		if(w1 <= 10.0 || w2 <= 10.0) continue;
		if(s == 0) continue;
		if(t == gr.num_vertices() - 1) continue;

		vector<int> v;
		v.push_back(s - 1);
		v.push_back(k - 1);
		v.push_back(t - 1);
		
		hs.add_node_list(v, 1);
	}
	return 0;
}

int scallop3::init_super_edges()
{
	mev.clear();
	med.clear();
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		vector<int> v;
		int s = (*it1)->source();
		v.push_back(s);
		mev.insert(PEV(*it1, v));
		med.insert(PED(*it1, 0));
	}
	return 0;
}

int scallop3::init_inner_weights()
{
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		double w = gr.get_edge_weight(e);
		edge_info ei = gr.get_edge_info(e);
		ei.weight = w;
		gr.set_edge_info(e, ei);
	}
	return 0;
}

int scallop3::init_vertex_map()
{
	v2v.clear();
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		v2v.push_back(i);
	}
	return 0;
}

int scallop3::decompose_tree(undirected_graph &ug, const vector<int> &u2e)
{
	undirected_graph ug2(ug);
	while(true)
	{
		int x = -1, y = -1;
		edge_iterator it1, it2;
		for(tie(it1, it2) = ug2.edges(); it1 != it2; it1++)
		{
			int s = (*it1)->source();
			int t = (*it1)->target();
			assert(s != t);
			if(s > t) 
			{
				s = (*it1)->target();
				t = (*it1)->source();
			}
			if(ug2.degree(s) != 1 && ug2.degree(t) != 1) continue;

			x = s;
			y = t;
			ug2.remove_edge(*it1);
			break;
		}

		if(x == -1 || y == -1) break;

		int xx = u2e[x];
		int yy = u2e[y];
		int e = merge_adjacent_edges(xx, yy);
		hs.replace(xx, yy, e);
		if(ug.degree(x) == 1) hs.replace(xx, e);
		if(ug.degree(y) == 1) hs.replace(yy, e);
	}

	assert(ug2.num_edges() == 0);

	for(int i = 0; i < u2e.size(); i++)
	{
		int e = u2e[i];
		hs.remove(e);
	}
	return 0;
}

int scallop3::decompose_trivial_vertex(int x)
{
	vector<int> u2e;
	undirected_graph ug;
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.in_edges(x); it1 != it2; it1++)
	{
		int e = e2i[*it1];
		u2e.push_back(e);
		ug.add_vertex();
	}
	for(tie(it1, it2) = gr.out_edges(x); it1 != it2; it1++)
	{
		int e = e2i[*it1];
		u2e.push_back(e);
		ug.add_vertex();
	}

	for(int i = 0; i < gr.in_degree(x); i++)
	{
		for(int j = 0; j < gr.out_degree(x); j++)
		{
			ug.add_edge(i, j + gr.in_degree(x));
		}
	}

	balance_vertex(x);
	decompose_tree(ug, u2e);
	return 0;
}

int scallop3::greedy_decompose(int num)
{
	for(int i = 1; i < gr.num_vertices() - 1; i++) balance_vertex(i);
	for(int i = 1; i < gr.num_vertices() - 1; i++) balance_vertex(i);

	smoother sm(gr);
	sm.smooth();

	int cnt = 0;
	int n1 = paths.size();
	while(true)
	{
		if(num != -1 && cnt >= num) break;

		VE v;
		double w = gr.compute_maximum_path_w(v);

		if(w <= 0.0) break;
		if(w <= transcript_min_expression) break;

		int e = split_merge_path(v, w);
		collect_path(e);
		cnt++;
	}
	int n2 = paths.size();
	printf("greedy decomposing produces %d / %d paths\n", n2 - n1, n2);
	return 0;
}

int scallop3::split_merge_path(const VE &p, double wx)
{
	vector<int> v;
	for(int i = 0; i < p.size(); i++)
	{
		assert(p[i] != null_edge);
		assert(e2i.find(p[i]) != e2i.end());
		v.push_back(e2i[p[i]]);
	}
	return split_merge_path(v, wx);
}

int scallop3::split_merge_path(const vector<int> &p, double ww)
{
	if(p.size() == 0) return -1;
	int ee = split_edge(p[0], ww);
	for(int i = 1; i < p.size(); i++)
	{
		int x = split_edge(p[i], ww);
		ee = merge_adjacent_equal_edges(ee, x);
	}
	return ee;
}

int scallop3::merge_adjacent_equal_edges(int x, int y)
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

	vector<int> v = mev[xx];
	v.insert(v.end(), mev[yy].begin(), mev[yy].end());

	if(mev.find(p) != mev.end()) mev[p] = v;
	else mev.insert(PEV(p, v));

	double sum1 = gr.get_in_weights(xt);
	double sum2 = gr.get_out_weights(xt);

	//printf("vertex = %d, sum1 = %.2lf, sum2 = %.2lf\n", xt, sum1, sum2);

	assert(fabs(sum1 - sum2) <= SMIN);

	double sum = (sum1 + sum2) * 0.5;
	double r1 = gr.get_vertex_weight(xt) * (wx0 + wy0) * 0.5 / sum;
	double r2 = gr.get_vertex_weight(xt) - r1;
	gr.set_vertex_weight(xt, r2);

	//printf("sum = %.2lf, wx0 + wy0 = %.2lf, r1 = %.2lf, r2 = %.2lf\n", sum, (wx0 + wy0) * 0.5, r1, r2);
	double reads = med[xx] + med[yy] + r1 * gr.get_vertex_info(xt).length;

	//printf("set reads %.2lf to new edge\n", reads);

	if(med.find(p) != med.end()) med[p] = reads;
	else med.insert(PED(p, reads));

	assert(i2e[n] == p);
	assert(e2i.find(p) != e2i.end());
	assert(e2i[p] == n);
	assert(e2i[i2e[n]] == n);

	remove_edge(x);
	remove_edge(y);

	return n;
}

int scallop3::remove_edge(int e)
{
	edge_descriptor ee = i2e[e];
	assert(ee != null_edge);
	int s = ee->source();
	int t = ee->target();

	e2i.erase(ee);
	i2e[e] = null_edge;
	gr.remove_edge(ee);

	return 0;
}

int scallop3::merge_adjacent_edges(int x, int y)
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
	int xy = merge_adjacent_equal_edges(x1, y1);

	return xy;
}

int scallop3::split_edge(int ei, double w)
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

	gr.set_edge_weight(ee, ww - w);		// old edge
	gr.set_edge_info(ee, eif);			// old edge
	gr.set_edge_weight(p2, w);			// new edge
	gr.set_edge_info(p2, eif);			// new edge

	if(mev.find(p2) != mev.end()) mev[p2] = mev[ee];
	else mev.insert(PEV(p2, mev[ee]));

	assert(med.find(ee) != med.end());
	double reads = med[ee];
	double reads1 = w / ww * reads;
	double reads2 = reads - reads1;
	
	//printf("ei = %d, w = %.2lf, ww = %.2lf, reads = %.2lf, reads1 = %.2lf, reads2 = %.2lf\n", ei, w, ww, reads, reads1, reads2);

	assert(reads1 >= 0);
	assert(reads2 >= 0);

	//printf("set reads %.2lf to new edge, %.2lf to old edge\n", reads1, reads2);

	med[ee] = reads2;
	if(med.find(p2) != med.end()) med[p2] = reads1;
	else med.insert(PED(p2, reads1));

	int n = i2e.size();
	i2e.push_back(p2);
	e2i.insert(PEI(p2, n));

	return n;
}

bool scallop3::balance_vertex(undirected_graph &ug, const vector<int> & u2e)
{
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

	// routes weight variables
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
	if(f != GRB_OPTIMAL)
	{
		delete model;
		delete env;
		return false;
	}

	for(int i = 0; i < wvars.size(); i++)
	{
		double w = wvars[i].get(GRB_DoubleAttr_X);
		gr.set_edge_weight(i2e[u2e[i]], w);
	}

	delete model;
	delete env;
	return true;
}

int scallop3::balance_vertex(int v)
{
	if(gr.degree(v) <= 0) return 0;

	edge_iterator it1, it2;
	double w1 = 0, w2 = 0;
	for(tie(it1, it2) = gr.in_edges(v); it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		w1 += w;
	}
	for(tie(it1, it2) = gr.out_edges(v); it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		w2 += w;
	}

	assert(w1 >= SMIN);
	assert(w2 >= SMIN);

	double r1 = (w1 > w2) ? 1.0 : w2 / w1;
	double r2 = (w1 < w2) ? 1.0 : w1 / w2;

	for(tie(it1, it2) = gr.in_edges(v); it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		gr.set_edge_weight(*it1, w * r1);
	}
	for(tie(it1, it2) = gr.out_edges(v); it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		gr.set_edge_weight(*it1, w * r2);
	}

	return 0;
}

int scallop3::split_vertex(int x, const vector<int> &xe, const vector<int> &ye)
{
	assert(x != 0);
	assert(x != gr.num_vertices() - 1);
	if(xe.size() <= 0) return 0;
	if(ye.size() <= 0) return 0;

	double w1 = 0, w2 = 0;
	for(int i = 0; i < xe.size(); i++) w1 += gr.get_edge_weight(i2e[xe[i]]);
	for(int i = 0; i < ye.size(); i++) w2 += gr.get_edge_weight(i2e[ye[i]]);
	double r1 = w1 / gr.get_in_weights(x);
	double r2 = w2 / gr.get_out_weights(x);
	double ww1 = (r1 + r2) * 0.5 * gr.get_vertex_weight(x);
	double ww2 = gr.get_vertex_weight(x) - ww1;

	int n = gr.num_vertices();
	assert(v2v.size() == n);

	gr.add_vertex();
	gr.set_vertex_info(n, gr.get_vertex_info(n - 1));
	gr.set_vertex_info(n - 1, gr.get_vertex_info(x));
	gr.set_vertex_weight(n, gr.get_vertex_weight(n - 1));
	gr.set_vertex_weight(n - 1, ww1);
	gr.set_vertex_weight(x, ww2);

	v2v.push_back(v2v[n - 1]);
	v2v[n - 1] = v2v[x];

	edge_iterator it1, it2;
	VE ve;
	for(tie(it1, it2) = gr.in_edges(n - 1); it1 != it2; it1++) ve.push_back(*it1);
	for(int i = 0; i < ve.size(); i++)
	{
		edge_descriptor e = ve[i];
		int s = e->source(); 
		int t = e->target();
		assert(t == n - 1);
		gr.move_edge(e, s, n);
	}
	assert(gr.degree(n - 1) == 0);

	for(int i = 0; i < xe.size(); i++)
	{
		edge_descriptor e = i2e[xe[i]];
		assert(e != null_edge);
		int s = e->source();
		int t = e->target();
		assert(t == x);
		gr.move_edge(e, s, n - 1);
	}

	for(int i = 0; i < ye.size(); i++)
	{
		edge_descriptor e = i2e[ye[i]];
		assert(e != null_edge);
		int s = e->source();
		int t = e->target();
		assert(s == x);
		gr.move_edge(e, n - 1, t);
	}

	return 0;
}

vector<int> scallop3::topological_sort()
{
	vector<PI> v;
	for(int i = 0; i < v2v.size(); i++)
	{
		v.push_back(PI(v2v[i], i));
	}
	sort(v.begin(), v.end());

	vector<int> vv;
	for(int i = 0; i < v.size(); i++)
	{
		vv.push_back(v[i].second);
	}

	return vv;
}

int scallop3::collect_existing_st_paths()
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

int scallop3::collect_path(int e)
{
	assert(mev.find(i2e[e]) != mev.end());
	assert(med.find(i2e[e]) != med.end());

	//printf("reads for edge %d = %.2lf\n", e, med[i2e[e]]);

	vector<int> v0 = mev[i2e[e]];
	vector<int> v;
	for(int i = 0; i < v0.size(); i++) v.push_back(v2v[v0[i]]);

	sort(v.begin(), v.end());

	int n = v2v[gr.num_vertices() - 1];
	assert(v[0] == 0);
	assert(v[v.size() - 1] < n);
	v.push_back(n);

	path p;
	p.reads = med[i2e[e]];
	p.abd = gr.get_edge_weight(i2e[e]);
	p.v = v;
	paths.push_back(p);

	gr.remove_edge(i2e[e]);
	e2i.erase(i2e[e]);
	i2e[e] = null_edge;

	return 0;
}

double scallop3::compute_smallest_edge(int x, int &e)
{
	e = -1;
	edge_iterator it1, it2;
	double sum1 = 0;
	double sum2 = 0;
	double ratio = DBL_MAX;
	for(tie(it1, it2) = gr.in_edges(x); it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		sum1 += w;
	}
	for(tie(it1, it2) = gr.out_edges(x); it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		sum2 += w;
	}

	assert(sum1 >= SMIN);
	assert(sum2 >= SMIN);
	for(tie(it1, it2) = gr.in_edges(x); it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		double r = w * 1.0 / sum1;
		if(r >= ratio) continue;
		ratio = r;
		e = e2i[*it1];
	}
	for(tie(it1, it2) = gr.out_edges(x); it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		double r = w * 1.0 / sum2;
		if(r >= ratio) continue;
		ratio = r;
		e = e2i[*it1];
	}
	assert(e >= 0);
	return ratio;
}

int scallop3::stats()
{
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		vector<PI> p = hs.get_routes(i, gr, e2i);
		router rt(i, gr, e2i, i2e, p);
		rt.build();
		rt.stats();
	}
	return 0;
}

int scallop3::print()
{
	int n = 0;
	for(int i = 0; i < gr.num_vertices(); i++) 
	{
		if(gr.degree(i) >= 1) n++;
	}

	int p1 = gr.compute_num_paths();
	int p2 = gr.compute_decomp_paths();
	printf("statistics: %lu edges, %d vertices, total %d paths, %d required\n", gr.num_edges(), n, p1, p2);

	//hs.print();

	if(output_tex_files == true)
	{
		draw_splice_graph(name + "." + tostring(round) + ".tex");
		//nested_graph nt(gr);
		//nt.draw(name + "." + tostring(round) + ".nt.tex");
	}

	printf("finish round %d\n\n", round);

	round++;

	return 0;
}

int scallop3::draw_splice_graph(const string &file) 
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
		//sprintf(buf, "%d:%.1lf:%d:%.2lf:%c", i, w, l, d, b);
		sprintf(buf, "%d:%.0lf:%d", i, w, l);
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
		//sprintf(buf, "%d:%.1lf:%d:%c", i, w, l, b);
		sprintf(buf, "%d:%.0lf", i, w);
		mes.insert(PES(i2e[i], buf));
	}
	
	vector<int> tp = topological_sort();
	gr.draw(file, mis, mes, 4.5, tp);
	return 0;
}
