#ifndef __SGRAPH_H__
#define __SGRAPH_H__

#include <vector>
#include "bundle.h"
#include "dgraph.h"

using namespace std;

// class for splice graph
class sgraph: public bundle
{
public:
	sgraph(const bbase &bb);
	~sgraph();

public:
	dgraph gr;				// splice graph
	MEI e2b;				// edge-descriptor to bridge index

public:
	int solve();
	int build();
	int print();
};

#endif
