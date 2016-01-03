#ifndef __BUNDLE_H__
#define __BUNDLE_H__

#include <stdint.h>
#include <cstring>
#include <string>
#include <vector>

#include <boost/icl/interval_map.hpp>

#include "boundary.h"
#include "hit.h"
#include "sam.h"

using namespace std;
using namespace boost;

typedef pair<size_t, size_t> PT;
typedef icl::right_open_interval<int32_t> ROI;
typedef icl::interval_map<int32_t, int32_t, icl::partial_absorber, less, icl::inplace_plus, icl::inter_section, ROI> imap_t;

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
	imap_t imap;					// interval map
	vector<boundary> boundaries;	// all types of boundaries
	vector<PT> bridges;				// splice bridges

public:
	int solve();

	int clear();
	int print();

	// add hit to hits
	int add_hit(bam_hdr_t *h, bam1_t *b);

	// build interval map
	int build_interval_map();
	int count_overlap_reads(int32_t p);

	// infer boundaries
	int infer_splice_boundaries();
	int infer_left_boundaries();
	int infer_right_boundaries();
	int add_start_boundary();
	int add_end_boundary();

	// build regions
	int build_regions();

	// check whether hits are sorted
	int check_left_ascending();
	int check_right_ascending();
};

#endif
