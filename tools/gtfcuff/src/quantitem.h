#ifndef __QUANTITEM_H__
#define __QUANTITEM_H__

#include <string>

using namespace std;

class quantitem
{
public:
	quantitem(const string &s);
	bool operator<(const quantitem &qt) const;

public:
	int parse(const string &s);

public:
	string transcript_id;
	int length;
	double elength;
	double tpm;
	double numreads;
};

#endif
