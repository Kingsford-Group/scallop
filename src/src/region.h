#ifndef __REGION_H__
#define __REGION_H__

#include <stdint.h>
#include <vector>
#include "common.h"

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

	int32_t asc_pos;				// ascending position
	int32_t desc_pos;				// descending position

public:
	int print();

public:
	int32_t locate_ascending_position();
	int32_t locate_descending_position();
};

#endif
