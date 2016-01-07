#ifndef __BBASE_H__
#define __BBASE_H__

#include <stdint.h>
#include <cstring>
#include <string>
#include <vector>

#include "hit.h"
#include "sam.h"

using namespace std;

class bbase
{
public:
	bbase();
	virtual ~bbase();

protected:
	int32_t tid;					// chromosome ID
	string chrm;					// chromosome name
	int32_t lpos;					// the leftmost boundary on reference
	int32_t rpos;					// the rightmost boundary on reference
	vector<hit> hits;				// hits

public:
	int add_hit(bam_hdr_t *h, bam1_t *b);
	int clear();

	int32_t get_tid();
	int32_t get_rpos();
	size_t get_num_hits();
};

#endif
