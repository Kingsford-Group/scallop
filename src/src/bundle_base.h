#ifndef __BUNDLE_BASE_H__
#define __BUNDLE_BASE_H__

#include <stdint.h>
#include <cstring>
#include <string>
#include <vector>

#include "hit.h"

using namespace std;

class bundle_base
{
public:
	bundle_base();
	virtual ~bundle_base();

protected:
	int32_t tid;					// chromosome ID
	string chrm;					// chromosome name
	int32_t lpos;					// the leftmost boundary on reference
	int32_t rpos;					// the rightmost boundary on reference
	vector<hit> hits;				// hits
	char strand;					// strandness
	double ave_isize;				// average of all insert size (excluding mapped portion)

public:
	int add_hit(const hit &ht);
	int clear();

	int set_chrm(const string &s);
	int32_t get_tid();
	int32_t get_rpos();
	size_t get_num_hits();
};

#endif
