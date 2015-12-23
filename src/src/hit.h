#ifndef __HIT_H__
#define __HIT_H__

#include "sam.h"

// this class is a C++ closure of bam1_t
class hit
{
public:
	hit(bam1_t *_b);
	~hit();

private:
	bam1_t *b;
};

#endif
