#include "hit.h"

hit::hit(bam1_t *b)
	:bam1_core_t(b->core)
{
}

hit::~hit()
{
}
