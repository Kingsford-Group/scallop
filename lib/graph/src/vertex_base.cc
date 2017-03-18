/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

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

PEEI vertex_base::in_edges() const
{
	return PEEI(si.begin(), si.end());
}

PEEI vertex_base::out_edges() const
{
	return PEEI(so.begin(), so.end());
}

int vertex_base::print() const
{
	printf("in-edges = ( ");
	for(edge_iterator it = si.begin(); it != si.end(); it++)
	{
		printf("[%d, %d] ", (*it)->source(), (*it)->target());
	}
	printf("), out-edges = ( ");
	for(edge_iterator it = so.begin(); it != so.end(); it++)
	{
		printf("[%d, %d] ", (*it)->source(), (*it)->target());
	}
	printf(")\n");
	return 0;
}
