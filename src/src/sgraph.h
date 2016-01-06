#ifndef __SGRAPH_H__
#define __SGRAPH_H__

#include <vector>
#include "common.h"

using namespace std;

// class for splice graph
class sgraph
{
public:
	sgraph();
	~sgraph();

public:
	dgraph gr;				// splice graph
	vector<double> ave;		// average weight for each vertex
	vector<double> dev;		// standard deviation for each vertex
	MED ewt;				// weights for edge

public:
	int add_node(double ave_abd, double dev_abd);
	int add_arc(int s, int t, double wt);
};

#endif
