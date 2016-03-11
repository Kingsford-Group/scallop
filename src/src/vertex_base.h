#ifndef __VERTEX_BASE_H__
#define __VERTEX_BASE_H__

#include <set>

using namespace std;

class edge_b;

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
	virtual int print() const;
};

#endif
