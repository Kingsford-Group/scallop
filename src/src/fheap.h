#ifndef __FHEAP_H__
#define __FHEAP_H__

#include "boost/heap/fibonacci_heap.hpp"

using namespace boost;
using namespace std;

typedef heap::fibonacci_heap<double, heap::compare<std::greater<double> > > fheap;
typedef fheap::handle_type handle_t;

int test_fheap();

#endif
