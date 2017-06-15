/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "scallop.h"
#include "config.h"

#include <cstdio>
#include <iostream>
#include <climits>
#include <cfloat>
#include <algorithm>

scallop::scallop()
{}

scallop::scallop(const string &s, const splice_graph &g, const hyper_set &h)
	: name(s), gr(g), hs(h)
{
	round = 0;
	if(output_tex_files == true) gr.draw(name + "." + tostring(round++) + ".tex");

	gr.get_edge_indices(i2e, e2i);
	//add_pseudo_hyper_edges();
	hs.build(gr, e2i);
	init_super_edges();
	init_vertex_map();
	init_inner_weights();
	init_nonzeroset();
}

scallop::~scallop()
{
}

int scallop::assemble()
{
	int c = classify();
	if(verbose >= 1) printf("process splice graph %s type = %d, vertices = %lu, edges = %lu\n", name.c_str(), c, gr.num_vertices(), gr.num_edges());

	//resolve_negligible_edges(false, max_decompose_error_ratio[NEGLIGIBLE_EDGE]);

	while(true)
	{	
		if(gr.num_vertices() > max_num_exons) break;

		bool b = false;

		b = resolve_trivial_vertex_fast(max_decompose_error_ratio[TRIVIAL_VERTEX]);
		if(b == true) continue;

		b = resolve_trivial_vertex(1, max_decompose_error_ratio[TRIVIAL_VERTEX]);
		if(b == true) continue;

		b = resolve_unsplittable_vertex(UNSPLITTABLE_SINGLE, 1, -0.5);
		if(b == true) continue;

		b = resolve_smallest_edges(max_decompose_error_ratio[SMALLEST_EDGE]);
		if(b == true) continue;

		b = resolve_negligible_edges(true, max_decompose_error_ratio[NEGLIGIBLE_EDGE]);
		if(b == true) continue;

		b = resolve_unsplittable_vertex(UNSPLITTABLE_MULTIPLE, 1, -0.5);
		if(b == true) continue;

		b = resolve_splittable_vertex(SPLITTABLE_HYPER, 1, max_decompose_error_ratio[SPLITTABLE_HYPER]);
		if(b == true) continue;

		b = resolve_unsplittable_vertex(UNSPLITTABLE_SINGLE, INT_MAX, -0.5);
		if(b == true) continue;

		b = resolve_unsplittable_vertex(UNSPLITTABLE_MULTIPLE, INT_MAX, -0.5);
		if(b == true) continue;

		b = resolve_unsplittable_vertex(UNSPLITTABLE_SINGLE, INT_MAX, max_decompose_error_ratio[UNSPLITTABLE_SINGLE]);
		if(b == true) continue;

		b = resolve_hyper_edge(2);
		if(b == true) continue;

		b = resolve_hyper_edge(1);
		if(b == true) continue;

		b = resolve_smallest_edges(DBL_MAX);
		if(b == true) continue;

		//summarize_vertices();

		b = resolve_trivial_vertex(2, max_decompose_error_ratio[TRIVIAL_VERTEX]);
		if(b == true) continue;

		break;
	}

	collect_existing_st_paths();
	greedy_decompose();

	//filter_transcripts();

	if(verbose >= 2) 
	{
		for(int i = 0; i < paths.size(); i++) paths[i].print(i);
		printf("finish assemble bundle %s\n\n", name.c_str());
	}

	return 0;
}

bool scallop::resolve_smallest_edges(double max_ratio)
{
	int se = -1;
	int root = -1;
	double ratio = max_ratio;
	bool flag = false;
	//for(set<int>::iterator it = nonzeroset.begin(); it != nonzeroset.end(); it++)
	//for(int i = 1; i < gr.num_vertices() - 1; i++)
	vector<int> vv(nonzeroset.begin(), nonzeroset.end());
	for(int k = 0; k < vv.size(); k++)
	{
		int i = vv[k];
		assert(gr.degree(i) >= 1);

		if(gr.in_degree(i) <= 1) continue;
		if(gr.out_degree(i) <= 1) continue;

		double r;
		int e = compute_smallest_edge(i, r);

		if(e == -1) continue;

		int s = i2e[e]->source();
		int t = i2e[e]->target();

		if(gr.out_degree(s) <= 1) continue;
		if(gr.in_degree(t) <= 1) continue;

		//if(hs.right_extend(e) || hs.left_extend(e)) continue; TODO
		if(hs.right_extend(e) && hs.left_extend(e)) continue;
		if(t == i && hs.right_extend(e)) continue;
		if(s == i && hs.left_extend(e)) continue;

		if(r < 0.01)
		{
			double w = gr.get_edge_weight(i2e[e]);
			if(verbose >= 2) printf("resolve small edge, edge = %d, weight = %.2lf, ratio = %.2lf, vertex = (%d, %d), degree = (%d, %d)\n", 
					e, w, r, s, t, gr.out_degree(s), gr.in_degree(t));

			remove_edge(e);
			hs.remove(e);
			flag = true;
			continue;
		}

		if(ratio < r) continue;

		ratio = r;
		se = e;
		root = i;
	}

	if(flag == true) return true;
	if(se == -1) return false;

	double sw = gr.get_edge_weight(i2e[se]);
	int s = i2e[se]->source();
	int t = i2e[se]->target();
	if(verbose >= 2) printf("resolve small edge, edge = %d, weight = %.2lf, ratio = %.2lf, vertex = (%d, %d), degree = (%d, %d)\n", 
			se, sw, ratio, s, t, gr.out_degree(s), gr.in_degree(t));

	remove_edge(se);
	hs.remove(se);

	return true;
}

bool scallop::resolve_negligible_edges(bool extend, double max_ratio)
{
	bool flag = false;
	//for(set<int>::iterator it = nonzeroset.begin(); it != nonzeroset.end(); it++)
	//for(int i = 1; i < gr.num_vertices() - 1; i++)
	vector<int> vv(nonzeroset.begin(), nonzeroset.end());
	for(int k = 0; k < vv.size(); k++)
	{
		int i = vv[k];
		assert(gr.in_degree(i) >= 1);
		assert(gr.out_degree(i) >= 1);
		if(gr.in_degree(i) <= 1) continue;
		if(gr.out_degree(i) <= 1) continue;

		double ww1 = gr.get_max_in_weight(i);
		double ww2 = gr.get_max_out_weight(i);

		set<int> s;
		edge_iterator it1, it2;
		for(tie(it1, it2) = gr.in_edges(i); it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			double w = gr.get_edge_weight(e);
			if(w > max_ratio * ww1) continue;
			if(extend && hs.right_extend(e2i[e])) continue;
			if(verbose >= 2) printf("resolve in-negligible edge, degree = (%d, %d), vertex = %d, weight = %.3lf / %.3lf\n", gr.in_degree(i), gr.out_degree(i), i, w, ww1);
			s.insert(e2i[e]);
		}
		for(tie(it1, it2) = gr.out_edges(i); it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			double w = gr.get_edge_weight(e);
			if(w > max_ratio * ww2) continue;
			if(extend && hs.left_extend(e2i[e])) continue;
			if(verbose >= 2) printf("resolve out-negligible edge, degree = (%d, %d), vertex = %d, weight = %.3lf / %.3lf\n", gr.in_degree(i), gr.out_degree(i), i, w, ww1);
			s.insert(e2i[e]);
		}

		for(set<int>::iterator it = s.begin(); it != s.end(); it++)
		{
			edge_descriptor e = i2e[*it];
			if(gr.out_degree(e->source()) <= 1) continue;
			if(gr.in_degree(e->target()) <= 1) continue;
			if(hs.right_extend(e2i[e]) && hs.left_extend(e2i[e])) continue;

			remove_edge(*it);
			hs.remove(*it);
			flag = true;
		}
	}

	return flag;
}

bool scallop::resolve_splittable_vertex(int type, int degree, double max_ratio)
{
	int root = -1;
	double ratio = max_ratio;
	vector<equation> eqns;
	//for(set<int>::iterator it = nonzeroset.begin(); it != nonzeroset.end(); it++)
	//for(int i = 1; i < gr.num_vertices() - 1; i++)
	vector<int> vv(nonzeroset.begin(), nonzeroset.end());
	for(int k = 0; k < vv.size(); k++)
	{
		int i = vv[k];
		assert(gr.degree(i) >= 1);
		if(gr.in_degree(i) <= 1) continue;
		if(gr.out_degree(i) <= 1) continue;

		assert(gr.in_degree(i) >= 1);
		assert(gr.out_degree(i) >= 1);

		MPII mpi = hs.get_routes(i, gr, e2i);
		router rt(i, gr, e2i, i2e, mpi);
		rt.classify();

		if(rt.type != type) continue;
		if(rt.degree > degree) continue;

		rt.build();
		assert(rt.eqns.size() == 2);

		//if(rt.degree == degree && ratio < rt.ratio) continue;
		if(ratio < rt.ratio) continue;

		root = i;
		ratio = rt.ratio;
		eqns = rt.eqns;
		//degree = rt.degree;
	}

	if(root == -1) return false;

	if(verbose >= 2) printf("resolve splittable vertex, type = %d, degree = %d, vertex = %d, ratio = %.2lf, degree = (%d, %d)\n", 
			type, degree, root, ratio, gr.in_degree(root), gr.out_degree(root));

	split_vertex(root, eqns[0].s, eqns[0].t);

	return true;
}

bool scallop::resolve_unsplittable_vertex(int type, int degree, double max_ratio)
{
	int root = -1;
	MPID pe2w;
	double ratio = max_ratio;
	bool flag = false;
	//for(int i = 1; i < gr.num_vertices() - 1; i++)
	//for(set<int>::iterator it = nonzeroset.begin(); it != nonzeroset.end(); it++)
	vector<int> vv(nonzeroset.begin(), nonzeroset.end());
	for(int k = 0; k < vv.size(); k++)
	{
		int i = vv[k];
		assert(gr.degree(i) >= 1);
		if(gr.in_degree(i) <= 1) continue;
		if(gr.out_degree(i) <= 1) continue;

		assert(gr.in_degree(i) >= 1);
		assert(gr.out_degree(i) >= 1);

		MPII mpi = hs.get_routes(i, gr, e2i);
		router rt(i, gr, e2i, i2e, mpi);
		rt.classify();

		if(rt.type != type) continue;
		if(rt.degree > degree) continue;

		rt.build();

		if(rt.ratio < -0.5)
		{
			if(verbose >= 2) printf("resolve unsplittable vertex, type = %d, degree = %d, vertex = %d, ratio = %.3lf, degree = (%d, %d)\n",
					type, degree, i, rt.ratio, gr.in_degree(i), gr.out_degree(i));
			decompose_vertex_extend(i, rt.pe2w);
			flag = true;
			continue;
		}

		if(rt.ratio > ratio) continue;

		root = i;
		ratio = rt.ratio;
		pe2w = rt.pe2w;
	}

	if(flag == true) return true;
	if(root == -1) return false;

	if(verbose >= 2) printf("resolve unsplittable vertex, type = %d, degree = %d, vertex = %d, ratio = %.3lf, degree = (%d, %d)\n",
			type, degree, root, ratio, gr.in_degree(root), gr.out_degree(root));

	decompose_vertex_extend(root, pe2w);
	return true;
}

bool scallop::resolve_hyper_edge(int fsize)
{
	edge_iterator it1, it2;
	vector<int> v1, v2;
	int root = -1;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		int e = e2i[*it1];
		int vs = (*it1)->source();
		int vt = (*it1)->target();

		MI s;
		s = hs.get_successors(e);
		//if(s.size() >= 2 && hs.right_extend(get_keys(s)) == false && (hs.left_extend(e) == false || gr.out_degree(vs) == 1))
		if(s.size() >= fsize && (hs.left_extend(e) == false || gr.out_degree(vs) == 1))
		{
			v1.push_back(e);
			v2 = get_keys(s);
			root = (*it1)->target();
			break;
		}

		s = hs.get_predecessors(e);
		//if(s.size() >= 2 && hs.left_extend(get_keys(s)) == false && (hs.right_extend(e) == false || gr.in_degree(vt) == 1))
		if(s.size() >= fsize && (hs.right_extend(e) == false || gr.in_degree(vt) == 1))
		{
			v1 = get_keys(s);
			v2.push_back(e);
			root = (*it1)->source();
			break;
		}
	}

	if(v1.size() == 0 || v2.size() == 0) return false;
	assert(v1.size() == 1 || v2.size() == 1);

	if(verbose >= 2) printf("resolve hyper edge, fsize = %d, vertex = %d, degree = (%d, %d), hyper edge = (%lu, %lu)\n",
			fsize, root, gr.in_degree(root), gr.out_degree(root), v1.size(), v2.size());

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

	double r1 = (sum1 < sum2) ? 1.0 : sum2 / sum1;
	double r2 = (sum1 > sum2) ? 1.0 : sum1 / sum2;

	for(int i = 0; i < w1.size(); i++) w1[i] *= r1;
	for(int i = 0; i < w2.size(); i++) w2[i] *= r2;

	set<int> ss;
	bool flag = false;
	for(int i = 0; i < w1.size(); i++)
	{
		for(int j = 0; j < w2.size(); j++)
		{
			double w = (w1[i] < w2[j]) ? w1[i] : w2[j];
			if(w < 1.0) continue;

			flag = true;
			int k1 = split_edge(v1[i], w);
			int k2 = split_edge(v2[j], w);
			int x = merge_adjacent_equal_edges(k1, k2);

			//printf(" split (%d, %d), w = %.2lf, weight = (%.2lf, %.2lf), merge (%d, %d) -> %d\n", v1[i], v2[j], w, w1[i], w2[j], k1, k2, x);

			hs.replace(v1[i], v2[j], x);
			if(k1 == v1[i]) hs.remove(v1[i]);
			if(k2 == v2[j]) hs.remove(v2[j]);
			//if(k1 == v1[i]) hs.replace(v1[i], x);
			//if(k2 == v2[j]) hs.replace(v2[j], x);
		}
	}
	return flag;
}

bool scallop::resolve_trivial_vertex(int type, double jump_ratio)
{
	int root = -1;
	double ratio = DBL_MAX;
	int se = -1;
	bool flag = false;
	//for(int i = 1; i < gr.num_vertices() - 1; i++)
	//for(set<int>::iterator it = nonzeroset.begin(); it != nonzeroset.end(); it++)
	vector<int> vv(nonzeroset.begin(), nonzeroset.end());
	for(int k = 0; k < vv.size(); k++)
	{
		int i = vv[k];
		assert(gr.degree(i) >= 1);
		if(gr.in_degree(i) <= 0) continue;
		if(gr.out_degree(i) <= 0) continue;
		assert(gr.in_degree(i) >= 1);
		assert(gr.out_degree(i) >= 1);

		if(gr.in_degree(i) >= 2 && gr.out_degree(i) >= 2) continue;
		if(classify_trivial_vertex(i, true) != type) continue;

		int e;
		double r = compute_balance_ratio(i);

		if(r < 1.02)
		{
			if(verbose >= 2) printf("resolve trivial vertex %d, type = %d, ratio = %.2lf, degree = (%d, %d)\n", i, type, 
					r, gr.in_degree(i), gr.out_degree(i));

			decompose_trivial_vertex(i);
			flag = true;
			continue;
		}

		if(ratio < r) continue;

		root = i;
		ratio = r;
		se = e;

		if(ratio < jump_ratio) break;
	}

	if(flag == true) return true;
	if(root == -1) return false;

	if(verbose >= 2) printf("resolve trivial vertex %d, type = %d, ratio = %.2lf, degree = (%d, %d)\n", root, type, 
			ratio, gr.in_degree(root), gr.out_degree(root));

	decompose_trivial_vertex(root);
	assert(gr.degree(root) == 0);
	return true;
}

bool scallop::resolve_trivial_vertex_fast(double jump_ratio)
{
	bool flag = false;
	//for(set<int>::iterator it = nonzeroset.begin(); it != nonzeroset.end(); it++)
	vector<int> vv(nonzeroset.begin(), nonzeroset.end());
	for(int k = 0; k < vv.size(); k++)
	{
		int i = vv[k];
		assert(gr.in_degree(i) >= 1);
		assert(gr.out_degree(i) >= 1);
		bool b = resolve_single_trivial_vertex_fast(i, jump_ratio);
		if(b == true) flag = true;
	}
	return flag;
}

bool scallop::resolve_single_trivial_vertex_fast(int i, double jump_ratio)
{
	if(gr.degree(i) == 0) return false;
	if(gr.in_degree(i) >= 2 && gr.out_degree(i) >= 2) return false;
	if(classify_trivial_vertex(i, false) != 1) return false;

	double r = compute_balance_ratio(i);
	if(r >= jump_ratio) return false;

	if(verbose >= 2) printf("resolve trivial vertex fast, vertex = %d, ratio = %.2lf, degree = (%d, %d)\n",
			i, r, gr.in_degree(i), gr.out_degree(i));

	decompose_trivial_vertex(i);
	assert(gr.degree(i) == 0);

	return true;
}

int scallop::summarize_vertices()
{
	for(set<int>::iterator it = nonzeroset.begin(); it != nonzeroset.end(); it++)
	{
		int i = (*it);
		assert(gr.in_degree(i) >= 1);
		assert(gr.out_degree(i) >= 1);
		
		if(gr.in_degree(i) == 1 || gr.out_degree(i) == 1)
		{
			int c = classify_trivial_vertex(i, false);
			printf("summary: trivial vertex %d, type = %d\n", i, c);
		}
		else
		{
			MPII mpi = hs.get_routes(i, gr, e2i);
			set<int> s1;
			set<int> s2;
			for(MPII::iterator it = mpi.begin(); it != mpi.end(); it++)
			{
				PI p = it->first;
				s1.insert(p.first);
				s2.insert(p.second);
			}
			router rt(i, gr, e2i, i2e, mpi);
			rt.classify();
			rt.build();
			printf("summary: nontrivial vertex %d, degree = (%d, %d), hyper edges = %lu, graph degree = (%lu, %lu), type = %d, degree = %d, ratio = %.3lf\n", 
					i, gr.in_degree(i), gr.out_degree(i), mpi.size(), s1.size(), s2.size(), rt.type, rt.degree, rt.ratio);
		}
	}
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

	//assert(p0 >= p1);
	bool b = (p0 <= p1) ? true : false;

	if(p0 == p1) return TRIVIAL;
	else return NORMAL;
}

int scallop::add_pseudo_hyper_edges()
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

int scallop::init_inner_weights()
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

int scallop::init_vertex_map()
{
	v2v.clear();
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		v2v.push_back(i);
	}
	return 0;
}

int scallop::init_nonzeroset()
{
	nonzeroset.clear();
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.degree(i) <= 0) continue;
		nonzeroset.insert(i);
	}
	return 0;
}

int scallop::decompose_vertex_extend(int root, MPID &pe2w)
{
	// compute degree of each edge
	map<int, int> mdegree;
	for(MPID::iterator it = pe2w.begin(); it != pe2w.end(); it++)
	{
		PI p = it->first;
		if(mdegree.find(p.first) == mdegree.end()) mdegree.insert(PI(p.first, 1));
		else mdegree[p.first]++;
		if(mdegree.find(p.second) == mdegree.end()) mdegree.insert(PI(p.second, 1));
		else mdegree[p.second]++;
	}

	// add edge-vertex for each adjacent edge of root
	// map adjacent edges of root to vertices [m, n)
	int m = gr.num_vertices() - 1;
	int n = m;
	map<int, int> ev1, ev2;
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.in_edges(root); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		int ei = e2i[e];
		int s = e->source();
		int t = e->target();
		assert(t == root);
		assert(mdegree.find(ei) != mdegree.end());
		if(mdegree[ei] >= 2) ev1.insert(PI(ei, n++));
	}
	for(tie(it1, it2) = gr.out_edges(root); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		int ei = e2i[e];
		int s = e->source();
		int t = e->target();
		assert(s == root);
		assert(mdegree.find(ei) != mdegree.end());
		if(mdegree[ei] >= 2) ev2.insert(PI(ei, n++));
	}

	// add vertices  and exchange sink
	for(int i = m; i < n; i++) 
	{
		gr.add_vertex();
		assert(nonzeroset.find(i) == nonzeroset.end());
		nonzeroset.insert(i);
		v2v.push_back(-1);
	}
	v2v[n] = v2v[m];
	exchange_sink(m, n);

	// set vertex info for new vertices
	// detach edge from root to new vertex
	for(map<int, int>::iterator it = ev1.begin(); it != ev1.end(); it++)
	{
		edge_descriptor e = i2e[it->first];
		int k = it->second;

		gr.move_edge(e, e->source(), k);

		double w = gr.get_edge_weight(e);
		gr.set_vertex_info(k, gr.get_vertex_info(root));
		gr.set_vertex_weight(k, w);
		v2v[k] = v2v[root];
	}
	for(map<int, int>::iterator it = ev2.begin(); it != ev2.end(); it++)
	{
		edge_descriptor e = i2e[it->first];
		int k = it->second;

		gr.move_edge(e, k, e->target());

		gr.set_vertex_info(k, vertex_info());
		gr.set_vertex_weight(k, 0);
		v2v[k] = -2;
	}

	// connecting edges according to pe2w
	for(MPID::iterator it = pe2w.begin(); it != pe2w.end(); it++)
	{
		int e1 = it->first.first;
		int e2 = it->first.second;
		double w = it->second;

		if(mdegree[e1] == 1)
		{
			assert(mdegree[e2] >= 2);
			assert(ev1.find(e1) == ev1.end());
			assert(ev2.find(e2) != ev2.end());
			edge_descriptor p = i2e[e1];
			int v1 = p->source();
			int v2 = ev2[e2];
			gr.move_edge(p, v1, v2);
			//gr.set_edge_weight(p, w);
		}
		else if(mdegree[e2] == 1)
		{
			assert(mdegree[e1] >= 2);
			assert(ev1.find(e1) != ev1.end());
			assert(ev2.find(e2) == ev2.end());
			edge_descriptor p = i2e[e2];
			int v1 = ev1[e1];
			int v2 = p->target();
			gr.move_edge(p, v1, v2);
			//gr.set_edge_weight(p, w);
		}
		else
		{
			assert(mdegree[e1] >= 2);
			assert(mdegree[e2] >= 2);

			assert(ev1.find(e1) != ev1.end());
			assert(ev2.find(e2) != ev2.end());
			int v1 = ev1[e1];
			int v2 = ev2[e2];

			edge_descriptor p = gr.add_edge(v1, v2);
		
			int z = i2e.size();
			i2e.push_back(p);
			e2i.insert(PEI(p, z));

			gr.set_edge_weight(p, w);
			gr.set_edge_info(p, edge_info());

			vector<int> v0;
			if(mev.find(p) != mev.end()) mev[p] = v0;
			else mev.insert(PEV(p, v0));

			hs.insert_between(e1, e2, z);
		}
	}

	assert(gr.degree(root) == 0);
	nonzeroset.erase(root);

	for(map<int, int>::iterator it = ev1.begin(); it != ev1.end(); it++)
	{
		int k = it->second;
		resolve_single_trivial_vertex_fast(k, max_decompose_error_ratio[TRIVIAL_VERTEX]);
	}
	for(map<int, int>::iterator it = ev2.begin(); it != ev2.end(); it++)
	{
		int k = it->second;
		resolve_single_trivial_vertex_fast(k, max_decompose_error_ratio[TRIVIAL_VERTEX]);
	}

	return 0;
}

int scallop::decompose_vertex_replace(int root, MPID &pe2w)
{
	// reassign weights
	MID md;
	for(MPID::iterator it = pe2w.begin(); it != pe2w.end(); it++)
	{
		int e1 = it->first.first;
		int e2 = it->first.second;
		double w = it->second;
		if(md.find(e1) == md.end()) md.insert(PID(e1, w));
		else md[e1] += w;
		if(md.find(e2) == md.end()) md.insert(PID(e2, w));
		else md[e2] += w;
	}
	for(MID::iterator it = md.begin(); it != md.end(); it++)
	{
		edge_descriptor e = i2e[it->first];
		double w = it->second;
		gr.set_edge_weight(e, w);
	}

	// assert that all edges are covered
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.in_edges(root); it1 != it2; it1++)
	{
		int e = e2i[*it1];
		assert(md.find(e) != md.end());
	}
	for(tie(it1, it2) = gr.out_edges(root); it1 != it2; it1++)
	{
		int e = e2i[*it1];
		assert(md.find(e) != md.end());
	}

	// remove hyper-edges that are not covered
	MPII mpi = hs.get_routes(root, gr, e2i);
	for(MPII::iterator it = mpi.begin(); it != mpi.end(); it++)
	{
		if(pe2w.find(it->first) != pe2w.end()) continue;
		hs.remove_pair(it->first.first, it->first.second);
	}

	map<int, int> m;
	for(MPID::iterator it = pe2w.begin(); it != pe2w.end(); it++)
	{
		int e1 = it->first.first;
		int e2 = it->first.second;
		if(m.find(e1) == m.end()) m.insert(PI(e1, 1));
		else m[e1]++;
		if(m.find(e2) == m.end()) m.insert(PI(e2, 1));
		else m[e2]++;
	}
	for(MPID::iterator it = pe2w.begin(); it != pe2w.end(); it++)
	{
		int e1 = it->first.first;
		int e2 = it->first.second;
		double w = it->second;

		int e = merge_adjacent_edges(e1, e2, w);

		hs.replace(e1, e2, e);

		if(m[e1] == 1) hs.replace(e1, e);
		if(m[e2] == 1) hs.replace(e2, e);
	}

	for(MPID::iterator it = pe2w.begin(); it != pe2w.end(); it++)
	{
		int e1 = it->first.first;
		int e2 = it->first.second;
		assert(hs.left_extend(e1) == false || hs.right_extend(e1) == false);
		assert(hs.left_extend(e2) == false || hs.right_extend(e2) == false);
		hs.remove(e1);
		hs.remove(e2);
	}

	assert(gr.degree(root) == 0);
	nonzeroset.erase(root);
	return 0;
}

int scallop::decompose_trivial_vertex(int x)
{
	balance_vertex(x);

	MPID pe2w;
	edge_iterator it1, it2;
	edge_iterator ot1, ot2;
	for(tie(it1, it2) = gr.in_edges(x); it1 != it2; it1++)
	{
		int e1 = e2i[*it1];
		double w1 = gr.get_edge_weight(*it1);
		for(tie(ot1, ot2) = gr.out_edges(x); ot1 != ot2; ot1++)
		{
			int e2 = e2i[*ot1];
			double w2 = gr.get_edge_weight(*ot1);
			double w = w1 <= w2 ? w1 : w2;

			pe2w.insert(PPID(PI(e1, e2), w));
		}
	}
	decompose_vertex_replace(x, pe2w);
	return 0;
}

int scallop::classify_trivial_vertex(int x, bool fast)
{
	int d1 = gr.in_degree(x);
	int d2 = gr.out_degree(x);
	if(d1 != 1 && d2 != 1) return -1;

	edge_iterator it1, it2;
	tie(it1, it2) = gr.in_edges(x);
	int e1 = e2i[*it1];
	tie(it1, it2) = gr.out_edges(x);
	int e2 = e2i[*it1];

	if(d1 == 1)
	{
		int s = i2e[e1]->source();
		if(gr.out_degree(s) == 1) return 1;
		if(fast && hs.right_dominate(e1) == true) return 1;
	}
	
	if(d2 == 1)
	{
		int t = i2e[e2]->target();
		if(gr.in_degree(t) == 1) return 1;
		if(fast && hs.left_dominate(e2) == true) return 1;
	}

	return 2;
}

int scallop::exchange_sink(int old_sink, int new_sink)
{
	VE ve;
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.in_edges(old_sink); it1 != it2; it1++) ve.push_back(*it1);

	for(int i = 0; i < ve.size(); i++)
	{
		edge_descriptor e = ve[i];
		int s = e->source(); 
		int t = e->target();
		assert(t == old_sink);
		gr.move_edge(e, s, new_sink);
	}
	assert(gr.degree(old_sink) == 0);
	return 0;
}

int scallop::split_merge_path(const VE &p, double wx)
{
	vector<int> v;
	for(int i = 0; i < p.size(); i++)
	{
		assert(p[i] != null_edge);
		assert(e2i.find(p[i]) != e2i.end());
		v.push_back(e2i[p[i]]);
		double w = gr.get_edge_weight(p[i]);
	}
	return split_merge_path(v, wx);
}

int scallop::split_merge_path(const vector<int> &p, double ww)
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
	assert(fabs(wx0 - wy0) <= SMIN);

	int lx1 = gr.get_edge_info(xx).length;
	int ly1 = gr.get_edge_info(yy).length;
	int lxt = gr.get_vertex_info(xt).length;
	int lxy = lx1 + ly1 + lxt;

	gr.set_edge_weight(p, wx0 * 0.5 + wy0 * 0.5);
	gr.set_edge_info(p, edge_info(lxy));

	vector<int> v = mev[xx];
	v.insert(v.end(), mev[yy].begin(), mev[yy].end());

	if(mev.find(p) != mev.end()) mev[p] = v;
	else mev.insert(PEV(p, v));

	double sum1 = gr.get_in_weights(xt);
	double sum2 = gr.get_out_weights(xt);
	double sum = (sum1 + sum2) * 0.5;
	double r1 = gr.get_vertex_weight(xt) * (wx0 + wy0) * 0.5 / sum;
	double r2 = gr.get_vertex_weight(xt) - r1;
	gr.set_vertex_weight(xt, r2);

	assert(i2e[n] == p);
	assert(e2i.find(p) != e2i.end());
	assert(e2i[p] == n);
	assert(e2i[i2e[n]] == n);

	remove_edge(x);
	remove_edge(y);

	if(gr.in_degree(xt) == 0 && gr.out_degree(xt) == 0)
	{
		assert(gr.degree(xt) == 0);
		nonzeroset.erase(xt);
	}

	return n;
}

int scallop::remove_edge(int e)
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

int scallop::merge_adjacent_edges(int x, int y, double ww)
{
	if(i2e[x] == null_edge) return -1;
	if(i2e[y] == null_edge) return -1;

	edge_descriptor xx = i2e[x];
	edge_descriptor yy = i2e[y];

	int xs = xx->source();
	int xt = xx->target();
	int ys = yy->source();
	int yt = yy->target();

	if(xt != ys) return merge_adjacent_edges(y, x, ww);
	assert(xt == ys);

	int x1 = split_edge(x, ww);
	int y1 = split_edge(y, ww);
	int xy = merge_adjacent_equal_edges(x1, y1);

	return xy;
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

	return merge_adjacent_edges(x, y, ww);
}

int scallop::split_edge(int ei, double w)
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

	int n = i2e.size();
	i2e.push_back(p2);
	e2i.insert(PEI(p2, n));

	return n;
}

int scallop::balance_vertex(int v)
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

	// use max-meature
	//double ww = (wv >= w1 && wv >= w2) ? wv : (w1 >= w2 ? w1 : w2);
	//assert(ww >= w1 && ww >= w2);

	// use sqrt-meature
	double ww = sqrt(w1 * w2);

	// use convex combination
	//double ww = sqrt(0.5 * w1 * w1 + 0.5 * w2 * w2);

	double r1 = ww / w1;
	double r2 = ww / w2;

	double m1 = 0, m2 = 0;
	for(tie(it1, it2) = gr.in_edges(v); it1 != it2; it1++)
	{
		double wx = gr.get_edge_weight(*it1);
		double wy = wx * r1;
		if(wy < 1.0)
		{
			m1 += (1.0 - wy);
			wy = 1.0;
		}
		gr.set_edge_weight(*it1, wy);
	}
	for(tie(it1, it2) = gr.out_edges(v); it1 != it2; it1++)
	{
		double wx = gr.get_edge_weight(*it1);
		double wy = wx * r2;
		if(wy < 1.0)
		{
			m2 += 1.0 - wy;
			wy = 1.0;
		}
		gr.set_edge_weight(*it1, wy);
	}

	if(m1 > m2)
	{
		edge_descriptor e = gr.max_out_edge(v);
		double w = gr.get_edge_weight(e);
		gr.set_edge_weight(e, w + m1 - m2);
	}
	else if(m1 < m2)
	{
		edge_descriptor e = gr.max_in_edge(v);
		double w = gr.get_edge_weight(e);
		gr.set_edge_weight(e, w + m2 - m1);
	}

	return 0;
}

double scallop::compute_balance_ratio(int v)
{
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

	//double ww = sqrt(0.5 * w1 * w1 + 0.5 * w2 * w2);
	double ww = 0.5 * w1 + 0.5 * w2;

	if(w1 >= w2) return w1 / w2;
	else return w2 / w1;
}

int scallop::split_vertex(int x, const vector<int> &xe, const vector<int> &ye)
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

	// vertex-n => new sink vertex
	// vertex-(n-1) => splitted vertex for xe and ye
	// vertex-x => splitted vertex for xe2 and ye2

	gr.add_vertex();
	assert(nonzeroset.find(n - 1) == nonzeroset.end());
	nonzeroset.insert(n - 1);
	gr.set_vertex_info(n, gr.get_vertex_info(n - 1));
	gr.set_vertex_info(n - 1, gr.get_vertex_info(x));
	gr.set_vertex_weight(n, gr.get_vertex_weight(n - 1));
	gr.set_vertex_weight(n - 1, ww1);
	gr.set_vertex_weight(x, ww2);

	v2v.push_back(v2v[n - 1]);
	v2v[n - 1] = v2v[x];

	// use vertex-n instead of vertex-(n-1) as sink vertex
	exchange_sink(n - 1, n);

	// attach edges in xe and ye to vertex-(n-1)
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

vector<int> scallop::topological_sort()
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
	assert(mev.find(i2e[e]) != mev.end());

	vector<int> v0 = mev[i2e[e]];
	vector<int> v;
	for(int i = 0; i < v0.size(); i++) 
	{
		if(v2v[v0[i]] < 0) continue;
		v.push_back(v2v[v0[i]]);
	}

	sort(v.begin(), v.end());

	int n = v2v[gr.num_vertices() - 1];
	assert(v[0] == 0);
	assert(v[v.size() - 1] < n);
	v.push_back(n);

	path p;
	p.abd = gr.get_edge_weight(i2e[e]);
	p.v = v;
	paths.push_back(p);

	gr.remove_edge(i2e[e]);
	e2i.erase(i2e[e]);
	i2e[e] = null_edge;

	return 0;
}

int scallop::greedy_decompose()
{
	if(gr.num_edges() == 0) return 0;

	for(int i = 1; i < gr.num_vertices() - 1; i++) balance_vertex(i);
	for(int i = 1; i < gr.num_vertices() - 1; i++) balance_vertex(i);

	int cnt = 0;
	int n1 = paths.size();
	while(true)
	{
		VE v;
		double w = gr.compute_maximum_path_w(v);
		if(w <= min_transcript_coverage) break;

		int e = split_merge_path(v, w);
		collect_path(e);
		cnt++;
	}
	int n2 = paths.size();
	if(verbose >= 2) printf("greedy decomposing produces %d / %d paths\n", n2 - n1, n2);
	return 0;
}

int scallop::filter_transcripts()
{
	double max_abd = 0;
	for(int i = 0; i < paths.size(); i++)
	{
		if(paths[i].abd < max_abd) continue;
		max_abd = paths[i].abd;
	}

	vector<path> v;
	for(int i = 0; i < paths.size(); i++)
	{
		if(paths[i].abd < min_transcript_coverage_ratio * max_abd) continue;
		v.push_back(paths[i]);
	}

	paths = v;
	return 0;
}

int scallop::compute_smallest_edge(int x, double &ratio)
{
	int e = -1;
	ratio = DBL_MAX;
	edge_iterator it1, it2;
	double sum1 = 0;
	double sum2 = 0;
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
	return e;
}

int scallop::print()
{
	/*
	int n0 = 0;
	int n1 = 0;
	for(int i = 1; i < gr.num_vertices() - 1; i++) 
	{
		if(gr.degree(i) == 0) n0++;
		if(gr.degree(i) >= 1) n1++;
	}
	assert(nonzeroset.size() == n1);
	*/

	printf("statistics: %lu edges, %lu vertices, %lu nonzero vertices\n", gr.num_edges(), gr.num_vertices(), nonzeroset.size());

	//int p1 = gr.compute_num_paths();
	//int p2 = gr.compute_decomp_paths();
	//printf("statistics: %lu edges, %d vertices, total %d paths, %d required\n", gr.num_edges(), n, p1, p2);

	//hs.print();

	//if(output_tex_files == true)
	//{
	//	draw_splice_graph(name + "." + tostring(round) + ".tex");
	//	nested_graph nt(gr);
	//	nt.draw(name + "." + tostring(round) + ".nt.tex");
	//}

	//printf("finish round %d\n\n", round);
	//round++;

	return 0;
}

int scallop::draw_splice_graph(const string &file) 
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
