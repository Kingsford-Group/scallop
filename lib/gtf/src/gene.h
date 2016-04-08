#ifndef __GTF_GENE_H__
#define __GTF_GENE_H__

#include <fstream>
#include <vector>

#include "exon.h"
#include "transcript.h"

using namespace std;

class gene
{
public:
	string seqname;
	string gene_id;
	vector<exon> exons;	
	vector<transcript> transcripts;			

public:
	int build_transcripts();
	int add_exon(const exon &ge);
	int write(ofstream &fout) const;	
	int print() const;
};

#endif
