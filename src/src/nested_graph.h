#ifndef __NESTED_GRAPH_H__
#define __NESTED_GRAPH_H__

#include "graph_base.h"

#include <map>
#include <cassert>

using namespace std;

class nested_graph : public graph_base
{
public:
	nested_graph();
	nested_graph(const graph_base &gr);
	virtual ~nested_graph();

public:
	virtual edge_descriptor add_edge(int s, int t);
	virtual int compute_out_ancestor(int v) const;
	virtual int compute_in_ancestor(int v) const;
};

#endif
