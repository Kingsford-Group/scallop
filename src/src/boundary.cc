#include <cstdio>
#include "boundary.h"

boundary::boundary()
{}

boundary::boundary(int _t, int32_t _p)
{
	type = _t;
	pos = _p;
	count = 0;
	min_qual = UINT32_MAX;
	max_qual = 0;
	score = 255;
}

boundary::boundary(int _t, int32_t _p, int32_t _c, uint32_t _min, uint32_t _max)
{
	type = _t;
	pos = _p;
	count = _c;
	min_qual = _min;
	max_qual = _max;
	score = 255;
}

boundary::boundary(const boundary &sp)
{
	type = sp.type;
	pos = sp.pos;
	count = sp.count;
	min_qual = sp.min_qual;
	max_qual = sp.max_qual;
	score = sp.score;
}

int boundary::print()
{
	printf("boundary: type = %1d, pos = %9d, count = %6d, min-qual = %4d, max-qual = %4d, score = %4d\n", type, pos, count, min_qual, max_qual, score);
	return 0;
}
