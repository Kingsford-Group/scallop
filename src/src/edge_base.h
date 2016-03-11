#ifndef __EDGE_BASE_H__
#define __EDGE_BASE_H__

#include <set>

using namespace std;

class edge_b
{
public:
	edge_b(int _s, int _t);

protected:
	int s;					// source
	int t;					// target

public:
	virtual int source() const;
	virtual int target() const;
	virtual int print() const;
};

typedef set<edge_b*>::iterator edge_iterator_b;
typedef pair<edge_iterator_b, edge_iterator_b> PEE_b;

#endif
