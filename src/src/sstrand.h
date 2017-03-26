/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __SSTRAND_H__
#define __SSTRAND_H__

#include "splice_graph.h"

using namespace std;

typedef pair<int, int> PI;

class sstrand
{
public:
	sstrand(const string &ch, splice_graph &gr);
	virtual ~sstrand();

public:
	string chrm;					// chrm name
	splice_graph &gr;				// splice graph
	vector<PI> segments;			// segments

public:
	int build();

private:
	int build_segments();
	int remove_inconsistent_edges();
	bool remove_inconsistent_edges(int k1, int k2);
	int analysis_segments();
	int analysis_segment(int k1, int k2);
	bool assign_vertex_strands();
	bool remove_inconsistent_strands();
	int analysis_strand();
	int refine_splice_graph();
};

#endif
