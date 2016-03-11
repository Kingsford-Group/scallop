#ifndef __PATH_SPACE_H__
#define __PATH_SPACE_H__

#include "boost_graph.h"
#include "algebra.h"

using namespace boost_graph;

class path_space
{
public:
	path_space();
	~path_space();

public:
	splice_graph gr;	// splice graph

public:
	VE i2e;
	MEI e2i;
	VVD pem0;			// backup
	VVD pem;			// path-edge-matrxi
	int rank;			// rank of pem

public:
	int solve(int n, int m);
	int solve(const string &file);
	int process();
	int simulate(int n, int m);
	int build_path_edge_matrix();
	int build_path_edge_matrix(int s, vector<int> &v);
	int sample_valid_paths();
	int sample_paths();
	int check_edges_cover();
	int check();

public:
	int print();
};

#endif
