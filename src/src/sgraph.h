#ifndef __SGRAPH_H__
#define __SGRAPH_H__

#include <vector>
#include "bundle.h"

using namespace std;

// class for splice graph
class sgraph: public bundle
{
public:
	sgraph();
	~sgraph();

public:
	dgraph gr;				// splice graph
	MEI e2b;				// edge-descriptor to bridge index

public:
	int solve();

	int build();
};

#endif
