#ifndef __GTF_GENE_H__
#define __GTF_GENE_H__

#include <fstream>
#include <vector>
#include <set>

#include "exon.h"
#include "transcript.h"

using namespace std;

class gene
{
public:
	vector<transcript> transcripts;			
	vector<exon> exons;

public:
	// build
	int add_exon(const exon &e);
	int build_transcripts();
	int add_transcript(const transcript &t);

	// modify
	int sort();
	int clear();
	int set_gene_id(const string &id);
	int remove_single_exon_transcripts();
	int remove_transcripts(double expression);
	int assign_RPKM(double factor);

	// fetch information
	string get_gene_id() const;
	string get_seqname() const;
	char get_strand() const;
	set<int32_t> get_exon_boundaries() const;
	PI32 get_bounds() const;

	// write
	int write(ofstream &fout) const;	
};

#endif
