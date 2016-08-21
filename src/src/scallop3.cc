#include "scallop3.h"
#include "config.h"
#include "smoother.h"

#include <cstdio>
#include <iostream>
#include <cfloat>

scallop3::scallop3()
{}

scallop3::scallop3(const string &s, const splice_graph &g)
	: name(s), gr(g)
{
	round = 0;
	if(output_tex_files == true) gr.draw(name + "." + tostring(round++) + ".tex");
	gr.get_edge_indices(i2e, e2i);
	init_super_edges();
	init_routers(gr.vhe);
	print();
}

scallop3::~scallop3()
{}

int scallop3::assemble()
{
	classify();
	decompose();
	return 0;
}

bool scallop3::decompose()
{
	int root = -1;
	double ratio = -1;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		router &rt = routers[i];
		double r = rt.divide();
		if(r < 0) continue;
		if(ratio >= 0 && ratio < r) continue;
		root = i;
		ratio = r;
	}

	if(root == -1) return false;

	printf("best ratio = %.2lf, root = %d\n", ratio, root);

	bool b = smooth_vertex(root, routers[root].eqns);
	if(b == false) return false;

	return 0;
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

int scallop3::init_routers(const vector<hyper_edge> &vhe)
{
	routers.clear();
	for(int i = 0; i < gr.num_vertices(); i++) routers.push_back(router(i, gr, e2i, i2e));

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

	if(ve.size() == 0)
	{
		printf("decompose trivial vertex %d\n", i);
	}

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

	routers[xs].split_out_edge(x, x1, ww / wx);
	routers[yt].split_in_edge(y, y1, ww / wy);

	int xy = merge_adjacent_equal_edges(x, y);

	routers[xs].replace_out_edge(x, xy);
	routers[yt].replace_in_edge(y, xy);
	routers[xt].remove_in_edge(x);
	routers[xt].remove_out_edge(y);

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

	return n;
}

bool scallop3::smooth_vertex(int v, const vector<equation> &eqns)
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
	gr.draw(file, mis, mes, 4.5);
	return 0;
}
