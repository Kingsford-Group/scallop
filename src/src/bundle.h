#ifndef __BUNDLE_H__
#define __BUNDLE_H__

#include <stdint.h>
#include <cstring>

#include "sam.h"

using namespace std;

class bundle
{
public:
	bundle();
	~bundle();

public:
	int32_t tid;				// chromosome ID
	char chr[1024];					// chromosome name
	int32_t lpos;				// the leftmost position
	int32_t rpos;				// the rightmost position
	bam1_t **hits;				// store hits
	int n_hits;					// number of hits
	int m_hits;					// size of the buffer

public:
	int add_hit(bam_hdr_t *h, bam1_t *b);
	int print();
};

#endif
