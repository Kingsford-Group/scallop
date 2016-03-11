#include "graph_base.h"

#include <cstdio>
#include <cassert>

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

int graph_b::add_edge(int s, int t)
{
	assert(s >= 0 && s < vv.size());
	assert(t >= 0 && t < vv.size());
	edge_b *e = new edge_b(s, t);
	assert(se.find(e) == se.end());
	se.insert(e);
	vv[s]->add_out_edge(e);
	vv[t]->add_in_edge(e);
	return 0;
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

PEE_b graph_b::edges()
{
	return PEE_b(se.begin(), se.end());
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

	PEE_b pb = gr.edges();
	for(edge_iterator_b it = pb.first; it != pb.second; it++)
	{
		gr.remove_edge(*it);
	}

	gr.print();

	return 0;
}
