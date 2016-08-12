#ifndef __REGION_H__
#define __REGION_H__

#include <stdint.h>
#include <vector>
#include "interval_map.h"
#include "partial_exon.h"

using namespace std;

typedef pair<int, int> PI;

class region
{
public:
	region(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype, const split_interval_map *_imap, const split_interval_map *_jmap);
	~region();

public:
	int32_t lpos;					// the leftmost boundary on reference
	int32_t rpos;					// the rightmost boundary on reference
	int ltype;						// type of the left boundary
	int rtype;						// type of the right boundary
	const split_interval_map *imap;	// pointer to match interval map
	const split_interval_map *jmap;	// pointer to indel interval map

	bool empty;						// whether this region is completely spliced
	int32_t lmid;					// >= 0, if reads starts from lpos
	int32_t rmid;					// >= 0, if reads ends at rpos

private:
	SIMI lit;
	SIMI rit;

public:
	int check_left_region();
	int check_right_region();
	bool empty_subregion(int32_t p1, int32_t p2);
	int print(int index) const;
	int build_partial_exons(vector<partial_exon> &pexons);
};

#endif
