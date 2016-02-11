#include "region.h"
#include "config.h"
#include "binomial.h"

region::region(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype, const imap_t *_imap)
	:lpos(_lpos), rpos(_rpos), imap(_imap), ltype(_ltype), rtype(_rtype)
{
	ave_abd = 0;
	dev_abd = 0;

	lit = imap->end();
	rit = imap->end();

	tie(lit, rit) = locate_boundary_iterators(*imap, lpos, rpos);

	check_empty();
	estimate_abundance();
}

region::region(const region &r)
	:lpos(r.lpos), rpos(r.rpos), imap(r.imap), ltype(r.ltype), rtype(r.rtype)
{
	ave_abd = r.ave_abd;
	dev_abd = r.dev_abd;
	empty = r.empty;
	lit = r.lit;
	rit = r.rit;
}

region& region::operator=(const region &r)
{
	lpos = r.lpos;
	rpos = r.rpos;
	imap = r.imap;
	ltype = r.ltype;
	rtype = r.rtype;
	lit = r.lit;
	rit = r.rit;
	return *this;
}

region::~region()
{}

int region::print(int index)
{
	char em = empty ? 'T' : 'F';
	char cl = left_break() ? 'T' : 'F';
	char cr = right_break() ? 'T' : 'F';
	printf("region %d: [%d-%d), type = (%d, %d), empty = %c, break = (%c, %c), ave-abd = %.1lf, std-abd = %.1lf\n",
			index, lpos, rpos, ltype, rtype, em, cl, cr, ave_abd, dev_abd);
	return 0;
}

int region::check_empty()
{
	if(ltype == RIGHT_SPLICE || rtype == LEFT_SPLICE)
	{
		empty = false;
		return 0;
	}

	if(lit == imap->end() || rit == imap->end())
	{
		empty = true;
		return 0;
	}

	int32_t s = compute_coverage(*imap, lit, rit);

	if(1.0 * s / (rpos - lpos) >= min_region_coverage) empty = false;
	else empty = true;

	return 0;
}


int region::estimate_abundance()
{
	if(empty == true) return 0;
	if(lit == imap->end() || rit == imap->end()) return 0;

	int sum = 0;
	for(ICI it = lit; ; it++)
	{
		assert(upper(it->first) > lower(it->first));
		sum += it->second * (upper(it->first) - lower(it->first));
		if(it == rit) break;
	}

	int len = upper(rit->first) - lower(lit->first);
	ave_abd = 1.0 * sum / len;

	double var = 0;
	for(ICI it = lit; ; it++)
	{
		assert(upper(it->first) > lower(it->first));
		var += (it->second - ave_abd) * (it->second - ave_abd) * (upper(it->first) - lower(it->first));
		if(it == rit) break;
	}

	dev_abd = sqrt(var / len);

	return 0;
}

bool region::left_break()
{
	if(ltype == RIGHT_SPLICE) return false;
	if(empty == true) return true;
	assert(lit != imap->end());
	if(lower(lit->first) > lpos) return true;
	else return false;
}

bool region::right_break()
{
	if(rtype == LEFT_SPLICE) return false;
	if(empty == true) return true;
	if(upper(rit->first) < rpos) return true;
	else false;
}
