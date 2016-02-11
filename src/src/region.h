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
	region(const region &r);
	region& operator=(const region &r);
	~region();

public:
	int32_t lpos;					// the leftmost boundary on reference
	int32_t rpos;					// the rightmost boundary on reference
	int ltype;						// type of the left boundary
	int rtype;						// type of the right boundary

	const imap_t *imap;				// pointer to a interval map
	ICI lit;						// left most interval, imap.end() if not valid
	ICI rit;						// right most interval, imap.end() if not valid

	bool empty;						// whether this region is completely spliced

	double ave_abd;					// average abundance
	double dev_abd;					// standard-deviation of abundance

public:
	int check_empty();
	int estimate_abundance();
	bool left_break();
	bool right_break();
	int print(int index);
};

#endif
