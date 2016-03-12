#ifndef __VERTEX_BASE_H__
#define __VERTEX_BASE_H__

#include <set>
#include "edge_base.h"

using namespace std;

class vertex_b
{
public:
	vertex_b();
	virtual ~vertex_b();

protected:
	set<edge_b*> si;		// in_edges
	set<edge_b*> so;		// out_edges

public:
	virtual int add_in_edge(edge_b *e);
	virtual int add_out_edge(edge_b *e);
	virtual int remove_in_edge(edge_b *e);
	virtual int remove_out_edge(edge_b *e);
	virtual int degree() const;
	virtual int in_degree() const;
	virtual int out_degree() const;
	virtual PEE_b in_edges() const;
	virtual PEE_b out_edges() const;
	virtual set<int> adjacent_vertices() const;
	virtual int print() const;
};

#endif
