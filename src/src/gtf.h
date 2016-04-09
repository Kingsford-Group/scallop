#ifndef __GTF_H__
#define __GTF_H__

#include <fstream>

#include "path.h"
#include "gene.h"
#include "interval_map.h"
#include "splice_graph.h"

using namespace std;

class gtf : public gene
{
public:
	gtf(const gene &g);

public:
	split_interval_map imap;

public:
	int build_splice_graph(splice_graph &gr);
	int output_gtf(ofstream &fout, const vector<path> &paths, const string &prefix) const;	
	int output_gtf(ofstream &fout) const;	

private:
	int build_split_interval_map();
	double compute_sum_expression();
	int add_vertices(splice_graph &gr);
	int add_edges(splice_graph &gr);
	int add_single_edge(int s, int t, double w, splice_graph &gr);
};

#endif
