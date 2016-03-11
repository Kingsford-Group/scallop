#ifndef __GTF_GENE_H__
#define __GTF_GENE_H__

#include <fstream>

#include "gtf_exon.h"
#include "interval_map.h"
#include "boost_graph.h"

using namespace std;
using namespace boost_graph;

class gtf_gene
{
public:
	vector<gtf_exon> exons;	
	vector< vector<int> > transcripts;			
	split_interval_map imap;

public:
	int build_splice_graph(splice_graph &gr);
	int output_gtf(ofstream &fout, const vector<path> &paths, const string &prefix) const;	
	int output_gtf(ofstream &fout) const;	
	int add_exon(const gtf_exon &ge);
	int print();

private:
	int build_transcripts();
	int build_split_interval_map();
	int32_t compute_sum_expression();
	int add_vertices(splice_graph &gr);
	int add_edges(splice_graph &gr);
	int add_single_edge(int s, int t, double w, splice_graph &gr);
};

#endif
