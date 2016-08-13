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
	bool operator< (const transcript &t) const;

public:
	string seqname;
	string source;
	string feature;
	string gene_id;

	string transcript_id;
	int32_t expression;
	vector<PI32> exons;
	char strand;

public:
	int add_exon(int s, int t);
	int add_exon(const exon &e);
	int sort();
	int write(ofstream &fout) const;
	int length() const;
	PI32 get_bounds() const;
	string label() const;
};

#endif
