#ifndef __REGION_H__
#define __REGION_H__

#include <stdint.h>
#include <vector>
#include "interval_map.h"

using namespace std;
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

	double ave_abd;					// average abundance
	double dev_abd;					// standard-deviation of abundance

private:
	vector<int> bins;			// average abundance for bins
	vector<int> s5end;			// score for 5end
	vector<int> s3end;			// score for 3end

public:
	int build();
	bool left_break() const;
	bool right_break() const;
	string label() const;
	int print(int index) const;
	int print_boundaries(int index) const;

private:
	int init_core_empty();
	int estimate_abundance();

	// identify boundaries inside this region
	int compute_bin_abundance();
	int compute_end_candidates();
};

#endif
