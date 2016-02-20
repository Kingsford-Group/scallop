#ifndef __GTF_GENE_H__
#define __GTF_GENE_H__

#include "gtf_exon.h"
#include "interval_map.h"
#include "splice_graph.h"

using namespace std;

class gtf_gene
{
public:
	vector<gtf_exon> exons;	
	vector< vector<int> > transcripts;			
	interval_map imap;

public:
	int build_splice_graph(splice_graph &gr);
	int add_exon(const gtf_exon &ge);
	int print();

private:
	int build_transcripts();
	int build_interval_map();
	int32_t compute_sum_expression();
	int add_vertices(splice_graph &gr);
	int add_edges(splice_graph &gr);
};

#endif
