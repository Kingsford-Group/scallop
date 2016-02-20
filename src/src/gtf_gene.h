#ifndef __GTF_GENE_H__
#define __GTF_GENE_H__

#include "gtf_exon.h"

using namespace std;

class gtf_gene
{
public:
	vector<gtf_exon> exons;	// all exons

public:
	int add_exon(const gtf_exon &ge);
};

#endif
