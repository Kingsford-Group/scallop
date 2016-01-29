#ifndef __SGRAPH_H__
#define __SGRAPH_H__

#include "bundle.h"
#include "dgraph.h"

// class for splice graph
class sgraph: public bundle
{
public:
	sgraph(const bbase &bb);
	virtual ~sgraph();

public:
	dgraph gr;						// splice graph

public:
	int solve();
	int build_graph();
	int check();
	int print(int index);
	int draw(const string &file);

	PEB get_max_in_edge(int x);
	PEB get_max_out_edge(int x);
};

#endif
