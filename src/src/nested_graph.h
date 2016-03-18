#ifndef __NESTED_GRAPH_H__
#define __NESTED_GRAPH_H__

#include "directed_graph.h"

#include <map>
#include <cassert>

using namespace std;

class nested_graph : public directed_graph
{
public:
	nested_graph();
	virtual ~nested_graph();

public:
	MEI mei;		// map the edge to the original vertex

public:
	int get_in_partner(int x); 
	int get_out_partner(int x);
	vector<int> get_pivots(const vector<int> &p);

public:
	int build(directed_graph &gr);
	int draw(const string &file);
};

#endif
