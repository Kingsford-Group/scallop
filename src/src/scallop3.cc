#include "scallop3.h"
#include "config.h"

#include <cstdio>
#include <iostream>
#include <cfloat>
#include <algorithm>

scallop3::scallop3()
{}

scallop3::scallop3(const string &s, const splice_graph &g)
	: name(s), gr(g)
{
	env = new GRBEnv();
	round = 0;
	if(output_tex_files == true) gr.draw(name + "." + tostring(round++) + ".tex");
	gr.get_edge_indices(i2e, e2i);
	init_super_edges();
	init_vertex_map();
	init_routers(gr.vhe);
	print();
}

scallop3::~scallop3()
{
	delete env;
}

int scallop3::assemble()
{
	classify();
	iterate();
	collect_existing_st_paths();
	return 0;
}

int scallop3::iterate()
{
	while(true)
	{
		bool b	= false;
		b = split_vertex();
		if(b == true) print();
		if(b == true) continue;

		b = decompose_trivial_vertex();
		if(b == true) print();
		if(b == true) continue;

		break;
	}
	return 0;
}

bool scallop3::split_vertex()
{
	int root = -1;
	double ratio = -1;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.in_degree(i) <= 1) continue;
		if(gr.out_degree(i) <= 1) continue;

		router &rt = routers[i];
		rt.update();
		rt.print();

		double r = rt.ratio;
		if(r < 0) continue;
		if(ratio >= 0 && ratio < r) continue;
		root = i;
		ratio = r;
	}

	if(root == -1) return false;
	if(ratio <= 0) return false;

	router &rt = routers[root];
	assert(rt.eqns.size() == 2);

	printf("split vertex %d, ratio = %.2lf, degree = (%d, %d), eqns = %lu\n", root, ratio, gr.in_degree(root), gr.out_degree(root), rt.eqns.size());
	for(int i = 0;i < rt.eqns.size(); i++) rt.eqns[i].print(99);

	split_vertex(root, rt.eqns[0].s, rt.eqns[0].t);

	return true;
}

bool scallop3::decompose_trivial_vertex()
{
	int root = -1;
	double ratio = -1;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.degree(i) == 0) continue;
		if(gr.in_degree(i) >= 2 && gr.out_degree(i) >= 2) continue;

		router &rt = routers[i];
		rt.update();

		double r = rt.ratio;
		if(r < 0) continue;
		if(ratio >= 0 && ratio < r) continue;
		root = i;
		ratio = r;
	}

	if(root == -1) return false;
	if(ratio <= 0) return false;

	router &rt = routers[root];
	assert(rt.eqns.size() == 1);

	printf("decompose trivial vertex %d, ratio = %.2lf, degree = (%d, %d), eqns = %lu\n", root, ratio, gr.in_degree(root), gr.out_degree(root), rt.eqns.size());
	for(int i = 0;i < rt.eqns.size(); i++) rt.eqns[i].print(99);

	balance_vertex(root);
	decompose_trivial_vertex(root);

	return true;
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

	assert(p0 >= p1);

	bool b = (p0 == p1) ? true : false;

	printf("\nprocess %s %s\n", name.c_str(), b ? "TRIVIAL" : "NORMAL");

	if(p0 == p1) return TRIVIAL;
	else return NORMAL;
}


int scallop3::init_super_edges()
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

int scallop3::init_vertex_map()
{
	v2v.clear();
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		v2v.push_back(i);
	}
	return 0;
}

int scallop3::init_routers(const vector<hyper_edge> &vhe)
{
	routers.clear();
	for(int i = 0; i < gr.num_vertices(); i++) routers.push_back(router(i, gr, e2i, i2e, env));

	for(int i = 0; i < vhe.size(); i++)
	{
		const hyper_edge &he = vhe[i];
		const vector<int> &vv = he.v;
		if(vv.size() <= 2) continue;

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
			int c = vhe[i].count;
			int x = vv[k + 1];

			if(c < min_router_count) continue;
			routers[x].add_route(PI(e1, e2), c);
		}
	}
	return 0;
}

bool scallop3::decompose_trivial_vertex(int i)
{
	vector<int> ve;
	return decompose_trivial_vertex(i, ve);
}

bool scallop3::decompose_trivial_vertex(int i, vector<int> &ve)
{
	if(i == 0) return false;
	if(i == gr.num_vertices() - 1) return false;
	if(gr.in_degree(i) >= 2 && gr.out_degree(i) >= 2) return false;
	if(gr.degree(i) == 0) return true;

	edge_iterator it1, it2, it3, it4;

	tie(it1, it2) = gr.in_edges(i);
	tie(it3, it4) = gr.out_edges(i);

	int e1 = e2i[*it1];
	int e2 = e2i[*it3];

	int d1 = gr.in_degree(i);
	int d2 = gr.out_degree(i);

	int ee = merge_adjacent_edges(e1, e2);

	//printf("degree = (%d, %d), e1 = %d, e2 = %d, ee = %d, degree' = (%d, %d)\n", d1, d2, e1, e2, ee,
	//		gr.in_degree(i), gr.out_degree(i));

	assert(gr.in_degree(i) < d1 || gr.out_degree(i) < d2);

	ve.push_back(ee);

	return decompose_trivial_vertex(i, ve);
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

	routers[xs].remove_out_edge(x);
	routers[xt].remove_in_edge(x);
	routers[ys].remove_out_edge(y);
	routers[yt].remove_in_edge(y);
	routers[xs].touched = true;
	routers[xt].touched = true;
	routers[yt].touched = true;

	return n;
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
	int xy = merge_adjacent_equal_edges(x, y);

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

	gr.set_edge_weight(ee, w);
	gr.set_edge_info(ee, eif);
	gr.set_edge_weight(p2, ww - w);
	gr.set_edge_info(p2, eif);

	if(mev.find(p2) != mev.end()) mev[p2] = mev[ee];
	else mev.insert(PEV(p2, mev[ee]));

	int n = i2e.size();
	i2e.push_back(p2);
	e2i.insert(PEI(p2, n));

	routers[s].remove_out_edge(ei);
	routers[t].remove_in_edge(ei);
	routers[s].touched = true;
	routers[t].touched = true;

	return n;
}

int scallop3::balance_vertex(int v)
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

	int n = gr.num_vertices();
	assert(v2v.size() == n);

	gr.add_vertex();
	gr.set_vertex_weight(n, gr.get_vertex_weight(n - 1));
	gr.set_vertex_info(n, gr.get_vertex_info(n - 1));
	gr.set_vertex_weight(n - 1, gr.get_vertex_weight(x));
	gr.set_vertex_info(n - 1, gr.get_vertex_info(x));

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

	routers.push_back(routers[n - 1]);
	routers[n - 1] = routers[x];
	routers[n - 1].root = n - 1;
	routers[n - 1].touched = true;
	routers[x].touched = true;

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
	vector<int> v0 = mev[i2e[e]];
	vector<int> v;
	for(int i = 0; i < v0.size(); i++) v.push_back(v2v[v0[i]]);

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

	//for(int i = 0; i < v2v.size(); i++) printf("%d->%d, ", i, v2v[i]);
	//printf("\n");

	printf("finish round %d\n\n", round);

	if(output_tex_files == true)
	{
		draw_splice_graph(name + "." + tostring(round) + ".tex");
		//nested_graph nt(gr);
		//nt.draw(name + "." + tostring(round) + ".nt.tex");
	}

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
