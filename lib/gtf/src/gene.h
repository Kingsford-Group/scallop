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
	vector<transcript> transcripts;			
	vector<exon> exons;

public:
	int add_exon(const exon &e);
	int build_transcripts();
	int add_transcript(const transcript &t);
	int clear();
	int sort();
	string get_gene_id() const;
	string get_seqname() const;
	int set_gene_id(const string &id);
	int write(ofstream &fout) const;	
};

#endif
