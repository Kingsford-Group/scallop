#ifndef __HIT_H__
#define __HIT_H__

#include "sam.h"

class hit: public bam1_core_t
{
public:
	hit(bam1_t *b);
	~hit();

public:
	uint32_t cigar[7];			// cigar, use samtools
	int32_t spos[7];			// splice positions
	int n_spos;					// number of splice positions

public:
	int print();
	int infer_splice_positions();
};

#endif
