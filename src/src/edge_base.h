#ifndef __EDGE_BASE_H__
#define __EDGE_BASE_H__

#include <set>

using namespace std;

#define null_edge NULL

class edge_base
{
public:
	edge_base(int _s, int _t);

protected:
	int s;					// source
	int t;					// target

public:
	virtual int move(int x, int y);
	virtual int source() const;
	virtual int target() const;
	virtual int print() const;
};

typedef edge_base* edge_descriptor;
typedef set<edge_base*>::iterator edge_iterator;
typedef pair<edge_descriptor, bool> PEB;
typedef pair<edge_iterator, edge_iterator> PEE;

#endif
