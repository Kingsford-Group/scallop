#include "nested_graph.h"

nested_graph::nested_graph()
{}

nested_graph::nested_graph(const graph_base &gr)
{
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		add_vertex();
	}

	edge_iterator it1, it2;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		add_edge(s, t);
	}
}

nested_graph::~nested_graph()
{}

edge_descriptor nested_graph::add_edge(int s, int t)
{
	edge_descriptor e = graph_base::add_edge(s, t);
	edge_iterator it1, it2;
	for(tie(it1, it2) = edges(); it1 != it2; it1++)
	{
		edge_descriptor ee = *it1;
		if(e == ee) continue;
		bool b = intersect(e, ee);
		if(b == false) continue;
		int ss = e->source() < ee->source() ? e->source() : ee->source();
		int tt = e->target() > ee->target() ? e->target() : ee->target();
		remove_edge(e);
		remove_edge(ee);
		add_edge(ss, tt);
		return null_edge;
	}
	return e;
}

int nested_graph::compute_out_ancestor(int x) const
{
	set<int> m;
	edge_iterator it1, it2;
	for(tie(it1, it2) = out_edges(x); it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		assert(s == x);
		if(m.find(t) == m.end()) m.insert(t);
	}

	for(set<int>::iterator it = m.begin(); it != m.end(); it++)
	{
		vector<int> v;
		bfs(*it, v);
		bool f = true;
		for(set<int>::iterator it2 = m.begin(); it2 != m.end(); it2++)
		{
			if(v[*it2] >= 1) 
			{
				f = false;
				break;
			}
		}
		if(f == true) return *it;
	}
	return -1;
}

int nested_graph::compute_in_ancestor(int x) const
{
	set<int> m;
	edge_iterator it1, it2;
	for(tie(it1, it2) = in_edges(x); it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		assert(t == x);
		if(m.find(s) == m.end()) m.insert(s);
	}

	for(set<int>::iterator it = m.begin(); it != m.end(); it++)
	{
		vector<int> v;
		bfs_reverse(*it, v);
		bool f = true;
		for(set<int>::iterator it2 = m.begin(); it2 != m.end(); it2++)
		{
			if(v[*it2] >= 1) 
			{
				f = false;
				break;
			}
		}
		if(f == true) return *it;
	}
	return -1;
}

