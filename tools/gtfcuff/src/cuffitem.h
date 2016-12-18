#ifndef __CUFFITEM_H__
#define __CUFFITEM_H__

#include <string>

using namespace std;

class cuffitem
{
public:
	cuffitem(const string &s);

public:
	int parse(const string &s);
	int print(int n, char c) const;

public:
	string ref_gene_id;
	string ref_transcript_id;
	string gene_id;
	string transcript_id;
	char code;
	int length;
	double coverage;
};

bool cuffitem_cmp_coverage(const cuffitem &x, const cuffitem &y);
bool cuffitem_cmp_length(const cuffitem &x, const cuffitem &y);

#endif
