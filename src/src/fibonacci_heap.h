#ifndef __FHEAP_H__
#define __FHEAP_H__

#include "boost/heap/fibonacci_heap.hpp"

using namespace boost;

struct fnode
{
	fnode(int _v, double _w) : v(_v), w(_w) {}
	bool operator<(const fnode & x) const { return w < x.w; }
	int v;
	double w;
};

typedef heap::fibonacci_heap<fnode> fibonacci_heap;
typedef fibonacci_heap::handle_type handle_t;

int test_fibonacci_heap();

#endif
