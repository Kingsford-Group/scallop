#ifndef __UTIL_H__
#define __UTIL_H__

#include <stdint.h>
#include <map>
#include <sstream>
#include <cassert>

using namespace std;

// macros: using int64_t for two int32_t
#define pack(x, y) (int64_t)((((int64_t)(x)) << 32) | ((int64_t)(y)))
#define high32(x) (int32_t)((x) >> 32)
#define low32(x) (int32_t)(((x) << 32) >> 32)

// definitions
typedef map<int32_t, int> MPI;
typedef pair<int32_t, int> PPI;
typedef pair<int, int> PI;

// common small functions
template<typename T>
string tostring(T t)
{
	ostringstream s;
	s << t;
	return s.str();
}

template<typename T>
T compute_overlap(const pair<T, T> &x, const pair<T, T> &y)
{
	assert(x.first <= x.second);
	assert(y.first <= y.second);
	if(x.first > y.first) return compute_overlap(y, x);
	assert(x.first <= y.first);
	if(y.first >= x.second) return x.second - y.first;
	if(x.second <= y.second) return x.second - y.first;
	else return y.second - y.first;
}

#endif