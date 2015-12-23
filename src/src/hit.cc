#include "hit.h"

hit::hit(bam1_t *_b)
{
	b = bam_dup1(_b);
}

hit::~hit()
{
	bam_destroy1(b);
}
