#ifndef __ITEM_H__
#define __ITEM_H__

#include <string>

using namespace std;

class item
{
public:
	item(const string &s);

public:
	int parse(const string &s);
	bool operator<(const item &ge) const;
	int print() const;

public:
	string gene_id;
	string transcript_id;
	char code;
	int length;
	double coverage;
};

#endif
