/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __QUANTIFIER_H__
#define __QUANTIFIER_H__

#include "hit.h"
#include "path.h"
#include "super_hit.h"
#include "directed_graph.h"

typedef pair<set<int>, super_hit> PSSH;
typedef map<set<int>, super_hit> MSSH;

class quantifier
{
public:
	quantifier(const vector<hit> &v, vector<path> &p);

public:
	const vector<hit> &hits;		// reference to hits
	vector<path> &paths;			// paths to be quantified
	MSSH super_hits;				// list of super_hits
	directed_graph pgr;				// path graph
	vector< vector<int> > pps;		// the paths for each vertex
	MEI e2i;						// edge map, from edge to index
	VE i2e;							// edge map, from index to edge

public:
	int quantify();
	int print();

private:
	int build_path_graph();
	int build_graph_index();
	int build_super_hits();
	int build_passing_paths();
	int assign_paths_to_super_hits();
	int assign_paths_to_super_hit(const set<int> &s, super_hit &h);
	bool check_valid_phasing(const set<int> &s);
};

#endif
