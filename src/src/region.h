#ifndef __REGION_H__
#define __REGION_H__

#include <stdint.h>
#include <vector>
#include "interval_map.h"
#include "slope.h"
#include "partial_exon.h"

using namespace std;

class region
{
public:
	region(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype, const split_interval_map *_imap);
	~region();

private:
	int32_t lpos;					// the leftmost boundary on reference
	int32_t rpos;					// the rightmost boundary on reference
	int ltype;						// type of the left boundary
	int rtype;						// type of the right boundary
	const split_interval_map *imap;	// pointer to a interval map

	int32_t lcore;					// left core position
	int32_t rcore;					// right core position
	bool empty;						// whether this region is completely spliced

	vector<int> bins;				// average abundance for bins
	vector<slope> slopes;			// slopes
	vector<partial_exon> pexons;	// partial exons;

public:
	int build();
	int print(int index) const;

private:
	int init();
	int estimate_abundance(int ll, int rr, double &ave, double &dev);
	double compute_deviation(const split_interval_map &sim);

	int compute_bin_abundances();
	int compute_slopes();
	int select_slopes();
	int refine_slope(slope &s);
	int build_partial_exons();
};

#endif
