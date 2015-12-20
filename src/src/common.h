#ifndef __COMMON_H__
#define __COMMON_H__

#include <stdint.h>

// pack two coordinates, each is represented by int32_t, into a int64_t
inline int64_t pack(int32_t x, int32_t y) 
{
	int64_t z = x;
	return ( (z << 32) | y);
}

// get the first/second 32bits of a int64_t
#define unpack1(x) ((int32_t)(x >> 32))
#define unpack2(x) ((int32_t)(x))

#endif
