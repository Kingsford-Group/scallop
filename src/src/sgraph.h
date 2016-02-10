#ifndef __SGRAPH_H__
#define __SGRAPH_H__

#include "bundle.h"
#include "dgraph.h"
#include "path.h"

// class for splice graph
class sgraph: public bundle
{
public:
	sgraph(const bbase &bb);
	virtual ~sgraph();

public:
	dgraph gr;						// splice graph

	vector<path> paths;				// transcripts

public:
	int solve();
	int build_graph();
	int check();
	int print(int index);
	int draw(const string &file);

	int build_paths();
	int compute_maximum_path(path &p);

	PEB get_max_in_edge(int x);
	PEB get_max_out_edge(int x);
};

#endif
