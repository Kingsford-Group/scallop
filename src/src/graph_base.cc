#include "graph_base.h"

#include <cstdio>
#include <cassert>
#include <tuple>

using namespace std;

graph_b::graph_b()
{}

graph_b::~graph_b()
{
	for(int i = 0; i < vv.size(); i++) delete vv[i];
	for(edge_iterator_b it = se.begin(); it != se.end(); it++)
	{
		delete (*it);
	}
}

int graph_b::add_vertex()
{
	vertex_b *v = new vertex_b;
	vv.push_back(v);
	return 0;
}

PEB_b graph_b::add_edge(int s, int t)
{
	assert(s >= 0 && s < vv.size());
	assert(t >= 0 && t < vv.size());
	edge_b *e = new edge_b(s, t);
	assert(se.find(e) == se.end());
	se.insert(e);
	vv[s]->add_out_edge(e);
	vv[t]->add_in_edge(e);
	return PEB_b(e, true);
}

PEB_b graph_b::edge(int s, int t) const
{
	assert(s >= 0 && s < vv.size());
	assert(t >= 0 && t < vv.size());
	PEE_b p = vv[s]->out_edges();
	for(edge_iterator_b it = p.first; it != p.second; it++)
	{
		int x = (*it)->target();
		if(x != t) continue;
		return PEB_b(*it, true);
	}
	return PEB_b(NULL, false);
}

int graph_b::remove_edge(edge_b *e)
{
	if(se.find(e) == se.end()) return -1;
	vv[e->source()]->remove_out_edge(e);
	vv[e->target()]->remove_in_edge(e);
	delete e;
	se.erase(e);
	return 0;
}

int graph_b::remove_edge(int s, int t)
{
	vector<edge_b*> v;
	PEE_b p = vv[s]->out_edges();
	for(edge_iterator_b it = p.first; it != p.second; it++)
	{
		if((*it)->target() != t) continue;
		v.push_back(*it);
	}
	for(int i = 0; i < v.size(); i++)
	{
		remove_edge(v[i]);
	}
	return 0;
}

int graph_b::clear_vertex(int x)
{
	vector<edge_b*> v;
	PEE_b pi = vv[x]->in_edges();
	PEE_b po = vv[x]->out_edges();
	for(edge_iterator_b it = pi.first; it != pi.second; it++)
	{
		v.push_back(*it);
	}
	for(edge_iterator_b it = po.first; it != po.second; it++)
	{
		v.push_back(*it);
	}
	for(int i = 0; i < v.size(); i++)
	{
		remove_edge(v[i]);
	}
	return 0;
}

int graph_b::clear()
{
	for(int i = 0; i < vv.size(); i++) delete vv[i];
	for(edge_iterator_b it = se.begin(); it != se.end(); it++)
	{
		delete (*it);
	}
	vv.clear();
	se.clear();
	return 0;
}

int graph_b::degree(int v) const
{
	return vv[v]->degree();
}

int graph_b::in_degree(int v) const
{
	return vv[v]->in_degree();
}

int graph_b::out_degree(int v) const
{
	return vv[v]->out_degree();
}

PEE_b graph_b::in_edges(int v) const
{
	return vv[v]->in_edges();
}

PEE_b graph_b::out_edges(int v) const
{
	return vv[v]->out_edges();
}

PEE_b graph_b::edges() const
{
	return PEE_b(se.begin(), se.end());
}


size_t graph_b::num_vertices() const
{
	return vv.size();
}

size_t graph_b::num_edges() const
{
	return se.size();
}

PAA_b graph_b::adjacent_vertices(int v) const
{
	return vv[v]->adjacent_vertices();
}

int graph_b::print() const
{
	printf("total %lu vertices, %lu edges\n", vv.size(), se.size());
	for(int i = 0; i < vv.size(); i++)
	{
		printf("vertex %d: ", i);
		vv[i]->print();
	}

	for(edge_iterator_b it = se.begin(); it != se.end(); it++)
	{
		(*it)->print();
	}
	return 0;
}

int graph_b::test()
{
	graph_b gr;
	gr.add_vertex();
	gr.add_vertex();
	gr.add_vertex();
	gr.add_vertex();
	gr.add_vertex();

	gr.add_edge(0, 1);
	gr.add_edge(0, 1);
	gr.add_edge(0, 1);
	gr.add_edge(1, 2);
	gr.add_edge(2, 3);
	gr.add_edge(2, 4);
	gr.add_edge(3, 4);

	gr.print();

	edge_iterator_b it1, it2;
	for(tie(it1, it2) = gr.edges(); it1 != it2; it1++)
	{
		gr.remove_edge(*it1);
	}

	gr.print();

	return 0;
}
