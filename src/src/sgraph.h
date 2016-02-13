#ifndef __SGRAPH_H__
#define __SGRAPH_H__

#include "bundle.h"
#include "dgraph.h"
#include "path.h"

// class for splice graph
class sgraph
{
public:
	sgraph();
	virtual ~sgraph();

public:
	dgraph gr;						// splice graph
	vector<path> paths;				// transcripts

public:
	int solve();
	int build(const bundle &bd);
	int draw(const string &file);
	int print();

private:
	int update_weights();

	int build_paths();
	int compute_maximum_forward_path(path &p);
	int compute_maximum_path(path &p);
	int subtract_path(const path &p);
};

#endif
