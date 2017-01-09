#ifndef __ESTIMATOR_H__
#define __ESTIMATOR_H__

#include "splice_graph.h"
#include "path.h"

class estimator
{
public:
	estimator(const splice_graph &g, const vector<path> &p);

public:
	splice_graph gr;					// splice graph
	vector<path> paths;					// predicted transcripts
	vector< set<int> > vs;				// paths go through each vertex

public:
	int estimate();
	double iterate();
	int build_path_indices();
	int build_path_lengths();
	int build_path_length(path &p);
	int print();
};

#endif
