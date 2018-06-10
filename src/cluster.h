/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __CLUSTER_H__
#define __CLUSTER_H__

#include "gene.h"
#include "undirected_graph.h"
#include "bundle_base.h"

typedef map< int, vector<int> > MIV;
typedef pair< int, vector<int> > PIV;

class cluster
{
public:
	cluster(const vector<transcript> &v);
	MIV miv;
	undirected_graph gr;

public:
	vector<transcript> trs;		// input transcripts
	vector<transcript> cct;		// transcripts w.r.t. clusters

public:
	int solve();
	int print();

private:
	int split_with_num_exons();
	int build_graph();
	int clustering();
	bool verify_equal(int x, int y);
	bool verify_subset(int x, int y);
};

#endif
