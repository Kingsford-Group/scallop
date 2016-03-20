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
	vector<PI> partners;

public:
	int get_in_partner(int x); 
	int get_out_partner(int x);

	// x \in [0, 2n): [0, n) -> into a vertex, [n, 2n) outof a vertex
	int bfs_search(int x, vector<int> &table, vector<int> &open); 
	bool link(int s, int t, vector<PI> &p);
	int test_linking();

	int build(directed_graph &gr);
	int draw(const string &file);

private:
	int init(directed_graph &gr);
	int build_nests(directed_graph &gr);
	int build_partners(directed_graph &gr);
};

#endif
