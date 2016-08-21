#include "scallop3.h"
#include "config.h"

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
	divide();
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

int scallop3::divide()
{
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		routers[i].divide();
	}
	return 0;
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
