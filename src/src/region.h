#ifndef __REGION_H__
#define __REGION_H__

#include <stdint.h>
#include <vector>

#include "hit.h"

using namespace std;

typedef vector<hit>::iterator HIT;

class region
{
public:
	region();
	~region();

public:
	int32_t lpos;			// the leftmost boundary on reference
	int32_t rpos;			// the rightmost boundary on reference

	HIT lit;				// left pointer to the hits vector
	HIT rit;				// right pointer to the hits vector
};

#endif
