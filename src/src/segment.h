#ifndef __SEGMENT_H__
#define __SEGMENT_H__

#include <stdint.h>
#include <vector>
#include "partial_exon.h"
#include "interval_map.h"

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
	vector<int> bins;				// average abundance for bins

public:
	int clear();
	int build();
	int add_partial_exon(int k, const partial_exon &pe, double w);
	int print(int index) const;
};

#endif
