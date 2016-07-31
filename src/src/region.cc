#include "region.h"
#include "config.h"
#include "util.h"
#include "binomial.h"
#include <algorithm>

using namespace std;

region::region(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype, const split_interval_map *_imap)
	:lpos(_lpos), rpos(_rpos), imap(_imap), ltype(_ltype), rtype(_rtype)
{
	lcore = lpos;
	rcore = rpos;

	init();
	build_bins();
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

int region::build_bins()
{
	bins.clear();
	if(empty == true) return 0;

	int32_t bsize = slope_bin_size;
	int bnum = ceil((rcore - lcore) * 1.0 / bsize);
	bins.resize(bnum, 0.0);

	SIMI lit, rit;
	tie(lit, rit) = locate_boundary_iterators(*imap, lpos, rpos);
	for(SIMI it = lit; ; it++)
	{
		int o = it->second;
		int32_t l = lower(it->first);
		int32_t r = upper(it->first);
		int bl = (l - lcore) / bsize;
		int br = (r - 1 - lcore) / bsize;

		assert(r > l);
		assert(br >= bl);
		assert(bl >= 0);
		assert(br < bins.size());

		for(int b = bl; b <= br; b++)
		{
			int x = lcore + b * bsize;
			int y = lcore + (b + 1) * bsize;
			if(b == bl) x = l;
			if(b == br) y = r;
			
			bins[b] += o * (y - x);
		}
		
		if(it == rit) break;
	}
	
	for(int i = 0; i < bins.size() - 1; i++)
	{
		bins[i] = bins[i] / bsize;
	}

	bins[bnum - 1] = bins[bnum - 1] / (rcore - bsize * bnum + bsize);

	return 0;
}

int region::print(int index) const
{
	int32_t lc = compute_overlap(*imap, lcore);
	int32_t rc = compute_overlap(*imap, rcore - 1);
	printf("region %d: empty = %c, pos = [%d, %d), core = [%d, %d), bins = %lu, coverage = (%d, %d)\n", 
			index, empty ? 'T' : 'F', lpos, rpos, lcore, rcore, bins.size(), lc, rc);
	return 0;
}
