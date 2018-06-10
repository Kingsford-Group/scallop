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

public:
	int quantify();
	int print();

private:
	int build_super_hits();
	int build_path_graph();
	bool check_valid_phasing(const set<int> &s);
};

#endif
