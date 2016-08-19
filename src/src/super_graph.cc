#include "super_graph.h"
#include <algorithm>

super_graph::super_graph()
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
		build_single_splice_graph(gr, s, index);
		subs.push_back(gr);
		index++;
	}
	return 0;
}

int super_graph::build_single_splice_graph(splice_graph &gr, const set<int> &ss, int index)
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

int super_graph::print()
{
	for(int i = 0; i < subs.size(); i++)
	{
		splice_graph &gr = subs[i];
		int32_t lpos = gr.get_vertex_info(0).lpos;
		int32_t rpos = gr.get_vertex_info(gr.num_vertices() - 1).rpos;

		printf("subgraph %d, #vertices = %lu / %lu, range = [%d, %d)\n", i, gr.num_vertices() - 2, root.num_vertices() - 2, lpos, rpos);
		// TODO
		//gr.print();
	}
	return 0;
}
