#include "vertex_base.h"
#include "edge_base.h"

#include <cstdio>
#include <cassert>
#include <cstdio>

using namespace std;

vertex_b::vertex_b()
{}

vertex_b::~vertex_b()
{}

int vertex_b::add_in_edge(edge_b *e)
{
	assert(si.find(e) == si.end());
	si.insert(e);
	return 0;
}

int vertex_b::add_out_edge(edge_b *e)
{
	assert(so.find(e) == so.end());
	so.insert(e);
	return 0;
}

int vertex_b::remove_in_edge(edge_b *e)
{
	assert(si.find(e) != si.end());
	si.erase(e);
	return 0;
}

int vertex_b::remove_out_edge(edge_b *e)
{
	assert(so.find(e) != so.end());
	so.erase(e);
	return 0;
}

int vertex_b::print() const
{
	printf("in-edges = ( ");
	for(edge_iterator_b it = si.begin(); it != si.end(); it++)
	{
		printf("%d ", (*it)->source());
	}
	printf("), out-edges = ( ");
	for(edge_iterator_b it = so.begin(); it != so.end(); it++)
	{
		printf("%d ", (*it)->target());
	}
	printf(")\n");
	return 0;
}
