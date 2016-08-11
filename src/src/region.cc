#include "region.h"
#include "config.h"
#include "util.h"
#include <algorithm>

using namespace std;

region::region(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype, const split_interval_map *_imap, const split_interval_map *_jmap)
	:lpos(_lpos), rpos(_rpos), imap(_imap), jmap(_jmap), ltype(_ltype), rtype(_rtype)
{
	lcore = lpos;
	rcore = rpos;

	init();
}

region::~region()
{}

int region::init()
{
	SIMI lit, rit;
	tie(lit, rit) = locate_boundary_iterators(*imap, lpos, rpos);

	if(lit == imap->end() || rit == imap->end())
	{
		empty = true;
		return 0;
	}

	if(ltype != RIGHT_SPLICE) lcore = lower(lit->first);
	if(rtype != LEFT_SPLICE) rcore = upper(rit->first);
	assert(rcore > lcore);

	if(lcore > lpos) ltype = START_BOUNDARY;
	if(rcore < rpos) rtype = END_BOUNDARY;

	if(ltype == RIGHT_SPLICE || rtype == LEFT_SPLICE) 
	{
		empty = false;
		return 0;
	}

	empty = true;
	int32_t m = compute_max_overlap(*imap, lit, rit);
	int32_t s = compute_sum_overlap(*imap, lit, rit);
	if(m > min_max_region_overlap) empty = false;
	if(1.0 * s / (rcore - lcore) > min_average_overlap) empty = false;

	return 0;
}

int region::print(int index) const
{
	int32_t lc = compute_overlap(*imap, lcore);
	int32_t rc = compute_overlap(*imap, rcore - 1);
	printf("region %d: empty = %c, type = (%d, %d), pos = [%d, %d), core = [%d, %d), boundary coverage = (%d, %d)\n", 
			index, empty ? 'T' : 'F', ltype, rtype, lpos, rpos, lcore, rcore, lc, rc);
	return 0;
}
