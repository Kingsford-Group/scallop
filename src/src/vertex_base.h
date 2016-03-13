#ifndef __VERTEX_BASE_H__
#define __VERTEX_BASE_H__

#include <set>
#include "edge_base.h"

using namespace std;

class vertex_base
{
public:
	vertex_base();
	virtual ~vertex_base();

protected:
	set<edge_base*> si;		// in_edges
	set<edge_base*> so;		// out_edges

public:
	virtual int add_in_edge(edge_base *e);
	virtual int add_out_edge(edge_base *e);
	virtual int remove_in_edge(edge_base *e);
	virtual int remove_out_edge(edge_base *e);
	virtual int degree() const;
	virtual int in_degree() const;
	virtual int out_degree() const;
	virtual PEE in_edges() const;
	virtual PEE out_edges() const;
	virtual set<int> adjacent_vertices() const;
	virtual int print() const;
};

#endif