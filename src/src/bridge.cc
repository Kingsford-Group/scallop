#include <cstdio>
#include "bridge.h"
#include "config.h"

bridge::bridge()
{}

bridge::bridge(int _t, int64_t _p)
{
	type = _t;
	lpos = high32(_p);
	rpos = low32(_p);
	count = 0;
	min_qual = UINT32_MAX;
	max_qual = 0;
	score = 255;
}

bridge::bridge(int _t, int64_t _p, int32_t _c, uint32_t _min, uint32_t _max)
{
	type = _t;
	lpos = high32(_p);
	rpos = low32(_p);
	count = _c;
	min_qual = _min;
	max_qual = _max;
	score = 255;
	lrgn = -1;
	rrgn = -1;
}

bridge::bridge(const bridge &sp)
{
	type = sp.type;
	lpos = sp.lpos;
	rpos = sp.rpos;
	count = sp.count;
	min_qual = sp.min_qual;
	max_qual = sp.max_qual;
	score = sp.score;

	lrgn = sp.lrgn;
	rrgn = sp.rrgn;
}

bool bridge::operator<(const bridge &x) const
{
	if(lpos <= x.lpos) return true;
	else return false;
}

int bridge::print()
{
	printf("bridge: type = %d, region = [%d, %d), %d -> %d, count = %d, min-qual = %d, max-qual = %d, score = %d\n", 
			type, lpos, rpos, lrgn, rrgn, count, min_qual, max_qual, score);
	return 0;
}
