#include "super_graph.h"
#include <algorithm>
#include <cfloat>

super_graph::super_graph(const splice_graph &gr, const hyper_set &hs)
	:root(gr), hyper(hs)
{}

super_graph::~super_graph()
{}

int super_graph::build()
{
	build_undirected_graph();
	split_splice_graph();
	return 0;
}

int super_graph::build_undirected_graph()
{
	ug.clear();
	for(int i = 0; i < root.num_vertices(); i++)
	{
		ug.add_vertex();
	}

	edge_iterator it1, it2;
	for(tie(it1, it2) = root.edges(); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		int s = e->source();
		int t = e->target();
		if(s == 0) continue;
		if(t == root.num_vertices() - 1) continue;

		ug.add_edge(s, t);
	}

	return 0;
}

int super_graph::split_splice_graph()
{
	vector< set<int> > vv = ug.compute_connected_components();
	int index = 0;
	for(int k = 0; k < vv.size(); k++)
	{
		set<int> &s = vv[k];
		if(s.size() == 1 && *(s.begin()) == 0) continue;
		if(s.size() == 1 && *(s.begin()) == root.num_vertices() - 1) continue;
		splice_graph gr;
		hyper_set hs;
		build_single_splice_graph(gr, hs, s, index);
		subs.push_back(gr);
		hss.push_back(hs);
		index++;
	}
	return 0;
}

int super_graph::build_single_splice_graph(splice_graph &gr, hyper_set &hs, const set<int> &ss, int index)
{
	//printf("build single splice graph with index = %d\n", index);
	gr.clear();
	vector<int> vv(ss.begin(), ss.end());
	sort(vv.begin(), vv.end());
	assert(vv.size() >= 1);

	int32_t lpos = root.get_vertex_info(vv[0]).lpos;
	int32_t rpos = root.get_vertex_info(vv[vv.size() - 1]).rpos;

	// vertices
	gr.add_vertex();
	vertex_info vi0;
	vi0.lpos = lpos;
	vi0.rpos = lpos;
	gr.set_vertex_weight(0, 0);
	gr.set_vertex_info(0, vi0);

	for(int i = 0; i < vv.size(); i++)
	{
		int k = vv[i];
		gr.add_vertex();
		gr.set_vertex_weight(i + 1, root.get_vertex_weight(k));
		gr.set_vertex_info(i + 1, root.get_vertex_info(k));
		PI p(index, i + 1);
		a2b.insert(pair<int, PI>(k, p));
		b2a.insert(pair<PI, int>(p, k));
	}

	gr.add_vertex();
	vertex_info vin;
	vin.lpos = rpos;
	vin.rpos = rpos;
	gr.set_vertex_weight(vv.size() + 1, 0);
	gr.set_vertex_info(vv.size() + 1, vin);

	// edges
	edge_iterator it1, it2;
	for(tie(it1, it2) = root.out_edges(0); it1 != it2; it1++)
	{
		int t = (*it1)->target();
		if(ss.find(t) == ss.end()) continue;
		assert(a2b.find(t) != a2b.end());
		assert(a2b[t].first == index);
		int y = a2b[t].second;

		edge_descriptor e = gr.add_edge(0, y);
		gr.set_edge_weight(e, root.get_edge_weight(*it1));
		gr.set_edge_info(e, root.get_edge_info(*it1));
	}

	int n = root.num_vertices() - 1;
	for(int i = 0; i < vv.size(); i++)
	{
		int s = vv[i];
		assert(s != 0 && s != n);
		assert(a2b.find(s) != a2b.end());
		assert(a2b[s].first == index);
		int x = a2b[s].second;

		for(tie(it1, it2) = root.out_edges(s); it1 != it2; it1++)
		{
			int t = (*it1)->target();
			assert(t == n || ss.find(t) != ss.end());
			assert(t == n || a2b.find(t) != a2b.end());
			assert(t == n || a2b[t].first == index);
			int y = ((t == n) ? gr.num_vertices() - 1 : a2b[t].second);

			edge_descriptor e = gr.add_edge(x, y);
			gr.set_edge_weight(e, root.get_edge_weight(*it1));
			gr.set_edge_info(e, root.get_edge_info(*it1));
		}
	}

	// hyper-set
	hs.clear();
	for(MVII::iterator it = hyper.nodes.begin(); it != hyper.nodes.end(); it++)
	{
		vector<int> v = it->first;
		int c = it->second;

		bool b = true;
		vector<int> vv;
		for(int k = 0; k < v.size(); k++)
		{
			if(ss.find(v[k]) == ss.end()) b = false;
			if(b == false) break;
			assert(a2b.find(v[k]) != a2b.end());
			assert(a2b[v[k]].first == index);
			int x = a2b[v[k]].second;
			vv.push_back(x);
		}

		if(b == false) continue;

		for(int i = 0; i < vv.size(); i++) vv[i]--;
		hs.add_node_list(vv, c);
	}

	return 0;
}

int super_graph::get_root_vertex(int s, int x) const
{
	assert(s >= 0 && s < subs.size());
	assert(x >= 0 && x < subs[s].num_vertices());
	if(x == 0) return 0;
	if(x == subs[s].num_vertices() - 1) return root.num_vertices() - 1;
	PI p(s, x);
	map<PI, int>::const_iterator it = b2a.find(p);
	assert(it != b2a.end());
	int v = it->second;
	assert(v > 0 && v < root.num_vertices() - 1);
	return v;
}

vector<int> super_graph::get_root_vertices(int s, const vector<int> &x) const
{
	vector<int> vv;
	for(int i = 0; i < x.size(); i++)
	{
		int v = get_root_vertex(s, x[i]);
		vv.push_back(v);
	}
	return vv;
}

int super_graph::build_maximum_path_graph(splice_graph &gr, undirected_graph &mg)
{
	mg.clear();
	int n = gr.num_vertices() - 1;
	if(gr.out_degree(0) <= 1) return 0;
	if(gr.in_degree(n) <= 1) return 0;

	edge_iterator it1, it2;
	vector<int> u2v;
	MI v2u;
	for(tie(it1, it2) = gr.out_edges(0); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		int s = e->target();
		assert(v2u.find(s) == v2u.end());
		v2u.insert(PI(s, v2u.size()));
		u2v.push_back(s);
		mg.add_vertex();
	}
	for(tie(it1, it2) = gr.in_edges(n); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		int t = e->source();
		assert(v2u.find(t) == v2u.end());
		v2u.insert(PI(t, v2u.size()));
		u2v.push_back(t);
		mg.add_vertex();
	}

	for(int i = 0; i < u2v.size(); i++)
	{
		int x = u2v[i];
		int y = -1;
		if(i < gr.out_degree(0)) 
		{
			compute_maximum_path1(gr, x, y);
			int j = v2u[y];
			assert(j >= gr.out_degree(0));
			assert(j < gr.out_degree(0) + gr.in_degree(n));
			mg.add_edge(i, j);
		}
		else if(i >= gr.out_degree(0) && i < gr.out_degree(0) + gr.in_degree(n))
		{
			compute_maximum_path2(gr, x, y);
			int j = v2u[y];
			assert(j >= 0);
			assert(j < gr.out_degree(0));
			mg.add_edge(j, i);
		}
		else assert(false);
	}

	return 0;
}

double super_graph::compute_maximum_path1(splice_graph &gr, int s, int &t)
{
	MED med;
	int n = gr.num_vertices() - 1;
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.in_edges(n); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		double w = gr.get_edge_weight(e);
		med.insert(PED(e, w));
		gr.set_edge_weight(e, DBL_MAX / 10.0);
	}

	VE ve;
	double ww = gr.compute_maximum_st_path_w(ve, s, n);
	assert(ve.size() >= 1);
	edge_descriptor e = ve[ve.size() - 1];
	assert(e->target() == n);
	t = e->source();

	for(MED::iterator it = med.begin(); it != med.end(); it++)
	{
		edge_descriptor e = it->first;
		double w = it->second;
		gr.set_edge_weight(e, w);
	}

	return ww;
}

double super_graph::compute_maximum_path2(splice_graph &gr, int t, int &s)
{
	MED med;
	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.out_edges(0); it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		double w = gr.get_edge_weight(e);
		med.insert(PED(e, w));
		gr.set_edge_weight(e, DBL_MAX / 10.0);
	}

	VE ve;
	double ww = gr.compute_maximum_st_path_w(ve, 0, t);
	assert(ve.size() >= 1);
	edge_descriptor e = ve[0];
	assert(e->source() == 0);
	s = e->target();

	for(MED::iterator it = med.begin(); it != med.end(); it++)
	{
		edge_descriptor e = it->first;
		double w = it->second;
		gr.set_edge_weight(e, w);
	}

	return ww;
}


int super_graph::print()
{
	int d0 = root.out_degree(0);
	int dn = root.in_degree(root.num_vertices() - 1);
	printf("super graph, %lu vertices, starting = %d, ending = %d, %lu subgraphs\n",
			root.num_vertices(), d0, dn, subs.size());
	for(int i = 0; i < subs.size(); i++)
	{
		splice_graph &gr = subs[i];
		d0 = gr.out_degree(0);
		dn = gr.in_degree(gr.num_vertices() - 1);

		int32_t lpos = gr.get_vertex_info(0).lpos;
		int32_t rpos = gr.get_vertex_info(gr.num_vertices() - 1).rpos;

		undirected_graph mg;
		build_maximum_path_graph(gr, mg);

		vector< set<int> > vv = mg.compute_connected_components();

		printf(" subgraph %d, #vertices = %lu / %lu, starting = %d, ending = %d, %lu components, range = [%d, %d)\n", 
				i, gr.num_vertices() - 2, root.num_vertices() - 2, d0, dn, vv.size(), lpos, rpos);
	}
	return 0;
}
