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
		printf("edge (%d, %d) and edge (%d, %d) intersect = %c\n", ee->source(), ee->target(), e->source(), e->target(), b ? 'T' : 'F');
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

int nested_graph::exchange(int x, int y, int z)
{
	assert(x >= 0 && x < num_vertices());
	assert(y >= 0 && y < num_vertices());
	assert(z >= 0 && z < num_vertices());
	assert(x == compute_in_ancestor(y));
	assert(z == compute_out_ancestor(y));

	edge_iterator it1, it2;
	set<edge_descriptor> m;
	for(tie(it1, it2) = out_edges(x); it1 != it2; it1++)
	{
		if(check_directed_path((*it1)->target(), y) == false) continue;
		m.insert(*it1);
	}
	for(tie(it1, it2) = in_edges(y); it1 != it2; it1++)
	{
		if(check_directed_path(x, (*it1)->source()) == false) continue;
		m.insert(*it1);
	}
	for(tie(it1, it2) = out_edges(y); it1 != it2; it1++)
	{
		if(check_directed_path((*it1)->target(), z) == false) continue;
		m.insert(*it1);
	}
	for(tie(it1, it2) = in_edges(z); it1 != it2; it1++)
	{
		if(check_directed_path(y, (*it1)->source()) == false) continue;
		m.insert(*it1);
	}

	for(set<edge_descriptor>::iterator it = m.begin(); it != m.end(); it++)
	{
		int s = (*it)->source();
		int t = (*it)->target();
		printf("s = %d, t = %d\n", s, t);
		if(s == x && t == y) move_edge(*it, y, z);
		else if(s == x) move_edge(*it, y, t);
		else if(t == y) move_edge(*it, s, z);
		else if(s == y && t == z) move_edge(*it, x, y);
		else if(s == y) move_edge(*it, x, t);
		else if(t == z) move_edge(*it, s, y);
		else assert(false);
	}
	return 0;
}
