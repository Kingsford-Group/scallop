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

bool boundary::operator<(const boundary &x) const
{
	if(pos <= x.pos) return true;
	else return false;
}

int boundary::print(int index) const
{
	printf("boundary %d: type = %d, pos = %d, count = %d, min-qual = %d, max-qual = %d, score = %d\n", 
			index, type, pos, count, min_qual, max_qual, score);
	return 0;
}
