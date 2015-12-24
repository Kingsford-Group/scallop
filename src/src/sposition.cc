#include <cstdio>
#include "sposition.h"

sposition::sposition()
{}

sposition::sposition(int32_t _p, int32_t _c, uint32_t _min, uint32_t _max)
{
	pos = _p;
	count = _c;
	min_qual = _min;
	max_qual = _max;
}

sposition::sposition(const sposition &sp)
{
	pos = sp.pos;
	count = sp.count;
	min_qual = sp.min_qual;
	max_qual = sp.max_qual;
}

int sposition::print()
{
	printf("Sposition: pos = %9d, count = %6d, min-qual = %3d, max-qual = %3d\n", pos, count, min_qual, max_qual);
	return 0;
}
