#ifndef __NESTED_GRAPH_H__
#define __NESTED_GRAPH_H__

#include "graph_base.h"

#include <map>
#include <cassert>

using namespace std;

class nested_graph : public graph_base
{
public:
	nested_graph(const graph_base &gr);
	virtual ~nested_graph();

public:
	MEI mei;		// map the edge to the original vertex

public:
	int build(const graph_base &gr);
	int draw(const string &file);
};

#endif
