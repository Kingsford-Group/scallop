#include <cstdio>
#include "position.h"

position::position(const position &p)
	:pos(p.pos)
{}

position::position(int32_t _p)
	:pos(_p)
{}

position::~position()
{}

splice_pos::splice_pos(int32_t _p)
	:position(_p)
{}

splice_pos::splice_pos(int32_t _p, int32_t _c, uint32_t _min, uint32_t _max)
	:position(_p)
{
	count = _c;
	min_qual = _min;
	max_qual = _max;
}

splice_pos::splice_pos(const splice_pos &sp)
	:position(sp)
{
	count = sp.count;
	min_qual = sp.min_qual;
	max_qual = sp.max_qual;
}

int splice_pos::print()
{
	printf("splice_pos: pos = %9d, count = %6d, min-qual = %3d, max-qual = %3d\n", pos, count, min_qual, max_qual);
	return 0;
}
