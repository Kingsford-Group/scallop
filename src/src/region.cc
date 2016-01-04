#include "region.h"
#include "config.h"

region::region(int32_t _lpos, int32_t _rpos, const imap_t *_imap)
	:lpos(_lpos), rpos(_rpos), imap(_imap)
{
	locate_ascending_position();
	locate_descending_position();
}

region::region(const region &r)
	:lpos(r.lpos), rpos(r.rpos), imap(r.imap)
{
	locate_ascending_position();
	locate_descending_position();
}

region& region::operator=(const region &r)
{
	lpos = r.lpos;
	rpos = r.rpos;
	imap = r.imap;
	asc_pos = r.asc_pos;
	desc_pos = r.desc_pos;
	return *this;
}

region::~region()
{}

int region::print()
{
	printf("region: [%d,%d), core = [%d, %d)\n", lpos, rpos, asc_pos, desc_pos);
	return 0;
}

int region::locate_ascending_position()
{
	asc_pos = lpos;

	int32_t x = lpos;
	int32_t y = x + ascending_step;
	int xy = cumulate_overlap(*imap, x, y);
	while(true)
	{
		int32_t z =  y + ascending_step;
		if(z > rpos) break;

		int yz = cumulate_overlap(*imap, y, z);

		uint32_t score = compute_binomial_score(xy + yz, 0.5, yz);
		if(score < min_ascending_score) break;
		
		asc_pos = y;

		x = y;
		y = z;
	}
	return 0;
}

int region::locate_descending_position()
{
	desc_pos = rpos;

	int32_t z = rpos;
	int32_t y = z - descending_step;
	int yz = cumulate_overlap(*imap, y, z);
	while(true)
	{
		int32_t x =  y - descending_step;
		if(x < lpos) break;

		int xy = cumulate_overlap(*imap, x, y);

		uint32_t score = compute_binomial_score(xy + yz, 0.5, xy);
		if(score < min_descending_score) break;
		
		desc_pos = y;

		z = y;
		y = x;
	}
	return 0;
}
