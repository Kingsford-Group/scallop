#ifndef __REGION_H__
#define __REGION_H__

#include <stdint.h>
#include <vector>
#include "interval_map.h"
#include "slope.h"
#include "partial_exon.h"

using namespace std;

typedef pair<int, int> PI;

class region
{
public:
	region(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype, const split_interval_map *_imap);
	~region();

public:
	int32_t lpos;					// the leftmost boundary on reference
	int32_t rpos;					// the rightmost boundary on reference
	int ltype;						// type of the left boundary
	int rtype;						// type of the right boundary
	const split_interval_map *imap;	// pointer to a interval map

	int32_t lcore;					// left core position
	int32_t rcore;					// right core position
	bool empty;						// whether this region is completely spliced

	vector<int> bins;				// average abundance for bins
	vector<bool> adjust;			// whether bin is adjusted
	vector<PI> boundaries;			// list of end5 boundaries (indices of bins)

public:
	int init();
	int build_bins();
	int add_boundary(int xi, int type);
	int print(int index) const;
	int build_partial_exons(vector<partial_exon> &pexons);
};

#endif
