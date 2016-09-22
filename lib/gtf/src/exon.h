#ifndef __GTF_EXON_H__
#define __GTF_EXON_H__

#include <string>
#include <stdint.h>

using namespace std;

class exon
{
public:
	exon(const string &s);
	exon(const string &_transcript_id, const string &_gene_id, int32_t _start, int32_t _end, int32_t _expression, double _coverage);

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
	double coverage;
};

#endif
