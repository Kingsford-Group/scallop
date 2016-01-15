#ifndef __REGION_H__
#define __REGION_H__

#include <stdint.h>
#include <vector>
#include "imap.h"

using namespace std;
class region
{
public:
	region(int32_t _lpos, int32_t _rpos, const imap_t *_imap);
	region(const region &r);
	region& operator=(const region &r);
	~region();

public:
	int32_t lpos;					// the leftmost boundary on reference
	int32_t rpos;					// the rightmost boundary on reference
	const imap_t *imap;				// pointer to a interval map

	bool empty;						// whether this region is completely spliced
	int32_t asc_pos;				// ascending position, inclusive
	int32_t desc_pos;				// descending position, exclusive

	double ave_abd;					// average abundance
	double dev_abd;					// standard-deviation of abundance

public:
	int print(int index);

public:
	int check_empty();
	int locate_ascending_position();
	int locate_descending_position();
	int estimate_abundance();
};

#endif
