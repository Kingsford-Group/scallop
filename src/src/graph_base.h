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
	virtual ~graph_b();

protected:
	vector<vertex_b*> vv;
	set<edge_b*> se;

public:
	virtual int add_vertex();
	virtual int add_edge(int s, int t);
	virtual int remove_edge(edge_b *e);
	virtual PEE_b edges();

public:
	virtual int print() const;
	static int test();
};

#endif
