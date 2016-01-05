#include "region.h"
#include "config.h"

region::region(int32_t _lpos, int32_t _rpos, const imap_t *_imap)
	:lpos(_lpos), rpos(_rpos), imap(_imap)
{
	asc_pos = lpos;
	desc_pos = rpos;
	check_empty();
	if(empty == false)
	{
		locate_ascending_position();
		locate_descending_position();
	}
}

region::region(const region &r)
	:lpos(r.lpos), rpos(r.rpos), imap(r.imap)
{
	asc_pos = r.asc_pos;
	desc_pos = r.desc_pos;
	empty = r.empty;
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
	char em = empty ? 'T' : 'F';
	char cl = lpos < asc_pos ? 'T' : 'F';
	char cr = rpos > desc_pos ? 'T' : 'F';
	printf("region: [%d-%d), empty = %c, check = (%c, %c), core = [%d-%d), origin-length = %d, core-length = %d\n",
			lpos, rpos, em, cl, cr, asc_pos, desc_pos, rpos - lpos, desc_pos - asc_pos);
	return 0;
}

int region::check_empty()
{
	int n = (rpos - lpos < num_sample_positions) ? (rpos - lpos) : num_sample_positions;
	int t = (rpos - lpos) / n;
	int s = cumulate_overlap(*imap, lpos, rpos, t);
	int r = (rpos - lpos) / t;
	double x = 1.0 * s / r;
	if(x < min_average_overlap) empty = true;
	else empty = false;
	return 0;
}

int region::locate_ascending_position()
{
	asc_pos = lpos;

	int32_t x = lpos;
	int32_t y = x + ascending_step;
	int xy = cumulate_overlap(*imap, x, y, 1);
	while(true)
	{
		int32_t z =  y + ascending_step;
		if(z > rpos) break;

		int yz = cumulate_overlap(*imap, y, z, 1);

		uint32_t score = compute_binomial_score(xy + yz, 0.5, yz);
		if(xy > 0 && score < min_ascending_score) break;
		
		asc_pos = y;

		x = y;
		y = z;
		xy = yz;
	}
	return 0;
}

int region::locate_descending_position()
{
	desc_pos = rpos;

	int32_t z = rpos;
	int32_t y = z - descending_step;
	int yz = cumulate_overlap(*imap, y, z, 1);

	//printf("descending "); print();
	while(true)
	{
		int32_t x =  y - descending_step;
		if(x < lpos) break;

		int xy = cumulate_overlap(*imap, x, y, 1);

		uint32_t score = compute_binomial_score(xy + yz, 0.5, xy);

		//printf("(%d,%d,%d) xy = %d, yz = %d, score = %d, min_descending_score = %d\n", x, y, z, xy, yz, score, min_descending_score);

		if(yz > 0 && score < min_descending_score) break;
		
		desc_pos = y;

		z = y;
		y = x;
		yz = xy;
	}
	//printf("\n");

	return 0;
}
