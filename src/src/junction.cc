#include <cstdio>
#include "junction.h"
#include "config.h"

junction::junction()
{}

junction::junction(int64_t _p)
{
	lpos = high32(_p);
	rpos = low32(_p);
	count = 0;
	min_qual = UINT32_MAX;
	max_qual = 0;
	score = 255;
}

junction::junction(int64_t _p, int32_t _c, uint32_t _min, uint32_t _max)
{
	lpos = high32(_p);
	rpos = low32(_p);
	count = _c;
	min_qual = _min;
	max_qual = _max;
	score = 255;
	lrgn = -1;
	rrgn = -1;
}

junction::junction(const junction &sp)
{
	lpos = sp.lpos;
	rpos = sp.rpos;
	count = sp.count;
	min_qual = sp.min_qual;
	max_qual = sp.max_qual;
	score = sp.score;

	lrgn = sp.lrgn;
	rrgn = sp.rrgn;
}

bool junction::operator<(const junction &x) const
{
	if(lpos <= x.lpos) return true;
	else return false;
}

int junction::print(int index) const
{
	printf("junction %d: region = [%d, %d), %d -> %d, count = %d, min-qual = %d, max-qual = %d, score = %d\n", 
			index, lpos, rpos, lrgn, rrgn, count, min_qual, max_qual, score);
	return 0;
}
