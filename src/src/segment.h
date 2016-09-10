#ifndef __SEGMENT_H__
#define __SEGMENT_H__

#include <stdint.h>
#include <vector>
#include "partial_exon.h"
#include "interval_map.h"
#include "slope.h"

using namespace std;

class segment
{
public:
	segment(const split_interval_map *_imap);
	~segment();

public:
	const split_interval_map *imap;	// pointer to a interval map
	vector<int> plist;				// partial exon index
	vector<partial_exon> pexons;	// exons in this segment
	vector<double> offsets;			// weight offset
	vector<int> nbins;				// #bins for each pexon
	vector<int> bsets;				// set of bins
	vector<slope> seeds;			// candidate slopes

public:
	int clear();
	int build();
	int add_partial_exon(int k, const partial_exon &pe, double w);
	int build_bin_sets();
	int build_bins(const partial_exon &pe, vector<int> &bins, int offset);
	int build_seeds();
	int print(int index) const;

private:
	int32_t get_left_position(int x);
	int32_t get_right_position(int x);
};

#endif
