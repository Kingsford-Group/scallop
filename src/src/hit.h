#ifndef __HIT_H__
#define __HIT_H__

#include <string>
#include <vector>

#include "config.h"
#include "sam.h"

using namespace std;

class hit: public bam1_core_t
{
public:
	hit(bam1_t *b);
	~hit();

public:
	int32_t rpos;							// right position mapped to reference [pos, rpos)
	string qname;							// query name
	uint32_t cigar[MAX_NUM_CIGAR];			// cigar, use samtools

public:
	int print();
	int get_splice_positions(vector<int64_t> &v);
	int get_matched_intervals(vector<int64_t> &v);
};

// for sorting
bool hit_compare_left(const hit &x, const hit &y);
bool hit_compare_right(const hit &x, const hit &y);

#endif
