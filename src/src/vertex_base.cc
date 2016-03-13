#include "vertex_base.h"
#include "edge_base.h"

#include <cstdio>
#include <cassert>
#include <cstdio>

using namespace std;

vertex_base::vertex_base()
{}

vertex_base::~vertex_base()
{}

int vertex_base::add_in_edge(edge_base *e)
{
	assert(si.find(e) == si.end());
	si.insert(e);
	return 0;
}

int vertex_base::add_out_edge(edge_base *e)
{
	assert(so.find(e) == so.end());
	so.insert(e);
	return 0;
}

int vertex_base::remove_in_edge(edge_base *e)
{
	assert(si.find(e) != si.end());
	si.erase(e);
	return 0;
}

int vertex_base::remove_out_edge(edge_base *e)
{
	assert(so.find(e) != so.end());
	so.erase(e);
	return 0;
}

int vertex_base::degree() const
{
	return in_degree() + out_degree();
}

int vertex_base::in_degree() const
{
	return (int)(si.size());
}

int vertex_base::out_degree() const
{
	return (int)(so.size());
}

PEE vertex_base::in_edges() const
{
	return PEE(si.begin(), si.end());
}

PEE vertex_base::out_edges() const
{
	return PEE(so.begin(), so.end());
}

set<int> vertex_base::adjacent_vertices() const
{
	set<int> x;
	/*
	for(edge_iterator it = si.begin(); it != si.end(); it++)
	{
		int s = (*it)->source();
		if(x.find(s) == x.end()) x.insert(s);
	}
	*/
	for(edge_iterator it = so.begin(); it != so.end(); it++)
	{
		int s = (*it)->target();
		if(x.find(s) == x.end()) x.insert(s);
	}
	return x;
}

int vertex_base::print() const
{
	printf("in-edges = ( ");
	for(edge_iterator it = si.begin(); it != si.end(); it++)
	{
		printf("%d ", (*it)->source());
	}
	printf("), out-edges = ( ");
	for(edge_iterator it = so.begin(); it != so.end(); it++)
	{
		printf("%d ", (*it)->target());
	}
	printf(")\n");
	return 0;
}
