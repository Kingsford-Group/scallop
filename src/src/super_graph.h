#ifndef __SUPER_GRAPH_H__
#define __SUPER_GRAPH_H__

#include "undirected_graph.h"
#include "splice_graph.h"

#include <map>
#include <cassert>

using namespace std;

class super_graph
{
public:
	super_graph();
	virtual ~super_graph();

public:
	int build();

public:
	splice_graph gr;
	undirected_graph ug;
	vector<splice_graph> subgraphs;

private:
	int build_undirected_graph();
	int split_splice_graph();
};

#endif
