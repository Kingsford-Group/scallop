#ifndef __GTF_EXON_H__
#define __GTF_EXON_H__

#include <string>
#include <stdint.h>

using namespace std;

class exon
{
public:
	exon(const string &s);

public:
	int parse(const string &s);
	bool operator<(const exon &ge) const;
	int print() const;
	int length() const;

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
