#ifndef __BUNDLE_H__
#define __BUNDLE_H__

#include <stdint.h>
#include <cstring>
#include <string>
#include <vector>

#include "boundary.h"
#include "hit.h"
#include "sam.h"

using namespace std;

class bundle
{
public:
	bundle();
	~bundle();

public:
	int32_t tid;					// chromosome ID
	string chrm;					// chromosome name
	int32_t lpos;					// the leftmost boundary on reference
	int32_t rpos;					// the rightmost boundary on reference
	vector<hit> hits;				// hits
	vector<boundary> boundaries;	// all types of boundaries

public:
	int solve();

	int clear();
	int print();

	int add_hit(bam_hdr_t *h, bam1_t *b);

	// check whether hits are sorted
	int check_left_ascending();
	int check_right_ascending();

	int count_prefix_left_hits();
	int count_suffix_left_hits();
	int count_prefix_right_hits();
	int count_suffix_right_hits();

	int infer_splice_boundaries();
	int infer_left_boundaries();
	int infer_right_boundaries();
	int add_start_boundary();
	int add_end_boundary();

	int build_regions();
};

#endif
