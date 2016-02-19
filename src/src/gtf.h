#ifndef __GTF_H__
#define __GTF_H__

#include <string>
#include <stdint.h>

using namespace std;

class gtf_line
{
public:
	gtf_line(const string &s);
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
