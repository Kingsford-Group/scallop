#ifndef __GTF_EXON_H__
#define __GTF_EXON_H__

#include "path.h"
#include <string>
#include <stdint.h>

using namespace std;

class gtf_exon
{
public:
	gtf_exon(const string &s);
	int print();

public:
	string seqname;
	string source;
	string feature;
	int32_t start;
	int32_t end;
	double score;
	char strand;
	int frame;
	string transcript_id;
	string gene_id;
};

#endif
