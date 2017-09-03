/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __GTF_ITEM_H__
#define __GTF_ITEM_H__

#include <string>
#include <stdint.h>

using namespace std;

class item
{
public:
	item(const string &s);

public:
	int parse(const string &s);
	bool operator<(const item &ge) const;
	int print() const;
	int length() const;

public:
	string seqname;
	string source;
	string feature;
	string gene_id;
	string transcript_id;
	string transcript_type;
	string gene_type;
	int32_t start;
	int32_t end;
	double score;
	char strand;
	int frame;
	double coverage;
	double FPKM;
	double RPKM;
	double TPM;
};

#endif
