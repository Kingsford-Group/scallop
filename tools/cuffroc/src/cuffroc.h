#ifndef __COMPARE_H__
#define __COMPARE_H__

#include "item.h"
#include <vector>

using namespace std;

class cuffroc
{
public:
	cuffroc(const string &file, int n);

public:
	vector<item> items;
	int refsize;

public:
	int read(const string &file);
	int solve();
};

#endif
