#ifndef __NESTED_GRAPH_H__
#define __NESTED_GRAPH_H__

#include "directed_graph.h"
#include "undirected_graph.h"

#include <map>
#include <cassert>

using namespace std;

typedef pair<int, int> PI;

class nested_graph : public directed_graph
{
public:
	nested_graph();
	virtual ~nested_graph();

public:
	vector<PI> partners;	// in/out partners for each vertex
	VE i2e;					// edge map
	MEI e2i;				// edge map
	vector<int>	parents;	// parent for each edge
	vector<bool> linkable;	// whether this edge can be moved to boundary

public:
	int build(directed_graph &gr);
	int draw(const string &file);
	int print();

private:
	int init(directed_graph &gr);
	int build_nests(directed_graph &gr);
	int build_partners(directed_graph &gr);
	int build_parents();
	int build_parents(int x, const vector<int> &order);
	int build_linkable();
	int build_linkable_forward(int x);
	int build_linkable_backward(int x);
	bool link(int xs, int xt, int ys, int yt, vector<PI> &p);
	int test_linking(directed_graph &gr);
};

#endif
