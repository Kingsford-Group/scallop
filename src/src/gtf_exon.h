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

public:
	int parse(const string &s);
	bool operator<(const gtf_exon &ge) const;
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
	int32_t expression;
};

#endif
