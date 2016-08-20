#include "super_graph.h"
#include <algorithm>

super_graph::super_graph()
{}

super_graph::~super_graph()
{}

int super_graph::build()
{
	build_hyper_edges();
	assign_edge_weights();
	identify_couple_edges();
	build_undirected_graph();
	split_splice_graph();
	return 0;
}

int super_graph::build_hyper_edges()
{
	for(int i = 0; i < vhe.size(); i++)
	{
		hyper_edge &he = vhe[i];
		vector<int> &vv = he.v;
		if(vv.size() <= 1) continue;

		he.ve.clear();
		for(int k = 0; k < vv.size() - 1; k++)
		{
			PEB p = root.edge(vv[k], vv[k + 1]);
			if(p.second == false) he.ve.push_back(null_edge);
			else he.ve.push_back(p.first);
		}
	}
	return 0;
}

int super_graph::assign_edge_weights()
{
	for(int i = 0; i < vhe.size(); i++)
	{
		VE &ve = vhe[i].ve;
		for(int k = 0; k < ve.size(); k++)
		{
			edge_descriptor &e = ve[k];
			if(e == null_edge) continue;
			edge_info ei = root.get_edge_info(e);
			ei.weight += vhe[i].count;
			root.set_edge_info(e, ei);
		}
	}

	/*
	edge_iterator it1, it2;
	for(tie(it1, it2) = root.edges(); it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		int32_t lpos = root.get_vertex_info(s).rpos;
		int32_t rpos = root.get_vertex_info(t).lpos;
		double w1 = root.get_edge_weight(*it1);
		double w2 = root.get_edge_info(*it1).weight;
		printf("edge (%d, %d), (%d, %d), weight1 = %.2lf weight2 = %.2lf\n", s, t, lpos, rpos, w1, w2);
	}
	*/
	return 0;
}

int super_graph::identify_couple_edges()
{
	typedef map< edge_descriptor, set<int> > MESI;
	typedef pair< edge_descriptor, set<int> > PESI;

	MESI m;
	for(int i = 0; i < vhe.size(); i++)
	{
		VE &ve = vhe[i].ve;
		for(int k = 0; k < ve.size(); k++)
		{
			edge_descriptor &e = ve[k];
			if(e == null_edge) continue;
			if(m.find(e) == m.end())
			{
				set<int> s;
				s.insert(i);
				m.insert(PESI(e, s));
			}
			else
			{
				m[e].insert(i);
			}
		}
	}

	for(int i = 1; i < root.num_vertices() - 1; i++)
	{
		if(root.in_degree(i) <= 1) continue;
		if(root.out_degree(i) <= 1) continue;

		edge_iterator it1, it2, it3, it4;
		for(tie(it1, it2) = root.in_edges(i); it1 != it2; it1++)
		{
			edge_descriptor e1 = (*it1);
			if(m.find(e1) == m.end()) continue;
			for(tie(it3, it4) = root.out_edges(i); it3 != it4; it3++)
			{
				edge_descriptor e2 = (*it3);
				if(m.find(e2) == m.end()) continue;

				vector<int> v1(m[e1].begin(), m[e1].end());
				vector<int> v2(m[e2].begin(), m[e2].end());
				sort(v1.begin(), v1.end());
				sort(v2.begin(), v2.end());

				/*
				printf("pair of edges: (");
				printv(v1);
				printf(") : (");
				printv(v2);
				printf(")\n");
				*/

				if(v1 != v2) continue;

				cps.push_back(PEE(e1, e2));
				//printf("couple edge (%d, %d) -> (%d, %d)\n", e1->source(), e1->target(), e2->source(), e2->target());
			}
		}
	}

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

	// hyper-edges
	for(int i = 0; i < vhe.size(); i++)
	{
		vector<int> &v = vhe[i].v;
		if(v.size() <= 2) continue;

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
		
		hyper_edge he(vv, vhe[i].count);
		gr.vhe.push_back(he);
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

	/*
	for(int i = 0; i < vhe.size(); i++)
	{
		vhe[i].print(i);
	}
	*/
	return 0;
}
