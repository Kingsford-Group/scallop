#include <cstring>
#include <cassert>
#include <cstdio>

#include "hit.h"

hit::hit(bam1_t *b)
	:bam1_core_t(b->core)
{
	// copy cigar
	assert(n_cigar <= 7);
	assert(n_cigar >= 1);
	memcpy(cigar, bam_get_cigar(b), 4 * n_cigar);
}

hit::~hit()
{
}

int hit::print()
{
	// print cigar
	if(n_cigar != 7) return 0;

	for(int i = 0; i < n_cigar; i++)
	{
		printf("cigar %d: op = %c, length = %3d\n", i, bam_cigar_opchr(cigar[i]), bam_cigar_oplen(cigar[i]));
	}
	return 0;
}

int hit::infer_splice_positions()
{
	return 0;
}
