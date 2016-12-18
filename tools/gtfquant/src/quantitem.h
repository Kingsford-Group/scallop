#ifndef __CUFFITEM_H__
#define __CUFFITEM_H__

#include <string>

using namespace std;

class quantitem
{
public:
	quantitem(const string &s);

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
