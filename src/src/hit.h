#ifndef __HIT_H__
#define __HIT_H__

#include <string>

#include "config.h"
#include "sam.h"

#define MAX_NUM_CIGAR 7
#define MIN_LEN_FLANK 1

using namespace std;

class hit: public bam1_core_t
{
public:
	hit(bam1_t *b);
	~hit();

public:
	int32_t rpos;							// right position mapped to reference [pos, rpos]
	string qname;							// query name
	uint32_t cigar[MAX_NUM_CIGAR];			// cigar, use samtools
	int32_t spos[MAX_NUM_CIGAR];			// splice positions
	int n_spos;								// number of splice positions

public:
	int print();
	int infer_splice_positions();
};

// for sorting
bool hit_compare_left(const hit &x, const hit &y);
bool hit_compare_right(const hit &x, const hit &y);

#endif
