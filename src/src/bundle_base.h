#ifndef __BUNDLE_BASE_H__
#define __BUNDLE_BASE_H__

#include <stdint.h>
#include <cstring>
#include <string>
#include <vector>

#include "hit.h"

using namespace std;

class bundle_base
{
public:
	bundle_base();
	virtual ~bundle_base();

protected:
	int32_t tid;					// chromosome ID
	string chrm;					// chromosome name
	int32_t lpos;					// the leftmost boundary on reference
	int32_t rpos;					// the rightmost boundary on reference
	vector<hit> hits;				// hits
	int phits;						// number of hits mapped to the positive strand
	int qhits;						// number of hits mapped to the reverse strand
	double ave_isize;				// average of all insert size (excluding mapped portion)

public:
	int add_hit(bam_hdr_t *h, bam1_t *b);
	int clear();

	int32_t get_tid();
	int32_t get_rpos();
	size_t get_num_hits();
};

#endif
