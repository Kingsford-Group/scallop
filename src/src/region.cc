#include "region.h"
#include "config.h"
#include "binomial.h"

region::region(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype, const interval_map *_imap)
	:lpos(_lpos), rpos(_rpos), imap(_imap), ltype(_ltype), rtype(_rtype)
{
	ave_abd = 0;
	dev_abd = 1;
	lcore = lpos;
	rcore = rpos;

	check_empty();
	estimate_abundance();
}

region::~region()
{}

int region::check_empty()
{
	SIMI lit, rit;
	tie(lit, rit) = locate_boundary_iterators(*imap, lpos, rpos);

	if(lit == imap->end() || rit == imap->end())
	{
		empty = true;
		return 0;
	}

	lcore = lower(lit->first);
	rcore = upper(rit->first);
	assert(rcore > lcore);

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

int region::estimate_abundance()
{
	if(empty == true) return 0;

	SIMI lit, rit;
	tie(lit, rit) = locate_boundary_iterators(*imap, lpos, rpos);

	ave_abd = 1.0 * compute_sum_overlap(*imap, lit, rit) / (rcore - lcore);

	double var = 0;
	for(SIMI it = lit; ; it++)
	{
		assert(upper(it->first) > lower(it->first));
		var += (it->second - ave_abd) * (it->second - ave_abd) * (upper(it->first) - lower(it->first));
		if(it == rit) break;
	}

	dev_abd = sqrt(var / (rcore - lcore));
	if(dev_abd < 0.1) dev_abd = 0.1;

	return 0;
}

bool region::left_break() const
{
	if(ltype == RIGHT_SPLICE) return false;
	if(empty == true) return true;
	if(lcore > lpos + 9) return true;
	else return false;
}

bool region::right_break() const
{
	if(rtype == LEFT_SPLICE) return false;
	if(empty == true) return true;
	if(rcore < rpos - 9) return true;
	else return false;
}

int region::print(int index) const
{
	char em = empty ? 'T' : 'F';
	char cl = left_break() ? 'T' : 'F';
	char cr = right_break() ? 'T' : 'F';
	printf("region %d: [%d-%d), type = (%d, %d), empty = %c, break = (%c, %c), core = [%d-%d), ave-abd = %.1lf, std-abd = %.1lf\n",
			index, lpos, rpos, ltype, rtype, em, cl, cr, lcore, rcore, ave_abd, dev_abd);
	return 0;
}
