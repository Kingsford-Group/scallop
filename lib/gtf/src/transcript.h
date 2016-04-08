#ifndef __GTF_TRANSCRIPT_H__
#define __GTF_TRANSCRIPT_H__

#include <fstream>
#include <string>
#include <vector>
#include "exon.h"

using namespace std;

typedef pair<int32_t, int32_t> PI32;

class transcript
{
public:
	transcript();
	~transcript();

public:
	vector<PI32> exons;
	string seqname;
	string source;
	string feature;
	string transcript_id;
	string gene_id;
	int32_t expression;

public:
	int add_exon(const exon &e);
	int sort();
	int write(ofstream &fout) const;
};

#endif
