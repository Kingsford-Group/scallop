#ifndef __REGION_H__
#define __REGION_H__

#include <stdint.h>
#include <vector>
#include "imap.h"

using namespace std;
class region
{
public:
	region(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype, const imap_t *_imap);
	~region();

public:
	int32_t lpos;					// the leftmost boundary on reference
	int32_t rpos;					// the rightmost boundary on reference
	int ltype;						// type of the left boundary
	int rtype;						// type of the right boundary

	const imap_t *imap;				// pointer to a interval map

	int32_t lcore;					// left core position
	int32_t rcore;					// right core position
	bool empty;						// whether this region is completely spliced

	double ave_abd;					// average abundance
	double dev_abd;					// standard-deviation of abundance

public:
	int estimate_abundance();
	bool left_break() const;
	bool right_break() const;
	int print(int index) const;
};

#endif
