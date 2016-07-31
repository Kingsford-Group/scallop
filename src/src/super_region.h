#ifndef __SUPER_REGION_H__
#define __SUPER_REGION_H__

#include <stdint.h>
#include <vector>
#include "region.h"
#include "slope.h"
#include "partial_exon.h"

using namespace std;

class super_region
{
public:
	super_region(const split_interval_map *_imap);
	~super_region();

private:
	const split_interval_map *imap;	// pointer to a interval map

	vector<region> regions;			// regions in this super-region
	vector<int> bins;				// average abundance for bins
	vector<slope> seeds;			// slopes candidates
	vector<slope> slopes;			// chosen slopes

public:
	vector<partial_exon> pexons;	// partial exons;

public:
	int clear();
	int build();
	int add_region(const region &r);
	int print(int index) const;
	int size() const;

private:
	int build_bins();
	int build_slopes();
	int build_slopes(int bin_num);
	int extend_slope(slope &s);
	int evaluate_slope(slope &s);
	int locate_bin(int x, int &xi, int &xb);
	int select_slopes(int si, int ti);
	int adjust_coverage(int si, int ti);
	int assign_boundaries();
	int build_partial_exons();
};

#endif
