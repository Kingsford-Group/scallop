#ifndef __HIT_H__
#define __HIT_H__

#include "sam.h"

#define MAX_NUM_CIGAR 7
#define MIN_LEN_FLANK 5


class hit: public bam1_core_t
{
public:
	hit(bam1_t *b);
	~hit();

public:
	uint32_t cigar[MAX_NUM_CIGAR];			// cigar, use samtools
	int32_t spos[MAX_NUM_CIGAR];			// splice positions
	int n_spos;								// number of splice positions

public:
	int print();
	int infer_splice_positions();
};

#endif
