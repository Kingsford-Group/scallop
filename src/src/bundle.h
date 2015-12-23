#ifndef __BUNDLE_H__
#define __BUNDLE_H__

#include <stdint.h>
#include <cstring>
#include <string>
#include <vector>

#include "hit.h"
#include "sam.h"

using namespace std;

class bundle
{
public:
	bundle();
	~bundle();

public:
	int32_t tid;				// chromosome ID
	string chr;					// chromosome name
	int32_t lpos;				// the leftmost position
	int32_t rpos;				// the rightmost position
	vector<hit> hits;			// hits

public:
	int add_hit(bam_hdr_t *h, bam1_t *b);
	int print();
	int clear();
};

#endif
