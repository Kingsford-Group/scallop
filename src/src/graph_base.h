#ifndef __GRAPH_BASE_H__
#define __GRAPH_BASE_H__

#include <vector>

#include "vertex_base.h"
#include "edge_base.h"

using namespace std;

class graph_b
{
public:
	graph_b();
	graph_b(const graph_b &gr);
	virtual ~graph_b();

protected:
	vector<vertex_b*> vv;
	set<edge_b*> se;

public:
	virtual int add_vertex();
	virtual PEB_b add_edge(int s, int t);
	virtual int remove_edge(edge_b *e);
	virtual int remove_edge(int s, int t);
	virtual int clear_vertex(int v);
	virtual int clear();

	virtual PEB_b edge(int s, int t) const;
	virtual size_t num_vertices() const;
	virtual size_t num_edges() const;
	virtual int degree(int v) const;
	virtual int in_degree(int v) const;
	virtual int out_degree(int v) const;
	virtual PEE_b edges() const;
	virtual PEE_b in_edges(int v) const;
	virtual PEE_b out_edges(int v) const;
	virtual set<int> adjacent_vertices(int v) const;

	virtual int move_edge(edge_b *e, int x, int y);

	virtual int print() const;
	static int test();
};

#endif
