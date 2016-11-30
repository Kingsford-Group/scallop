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
	bool operator<(const cuffitem &ge) const;
	int print(int n) const;

public:
	string ref_gene_id;
	string ref_transcript_id;
	string gene_id;
	string transcript_id;
	char code;
	int length;
	double coverage;
};

#endif
