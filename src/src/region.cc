#include "region.h"
#include "config.h"
#include "binomial.h"

region::region(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype, const imap_t *_imap)
	:lpos(_lpos), rpos(_rpos), imap(_imap), ltype(_ltype), rtype(_rtype)
{
	asc_pos = lpos;
	desc_pos = rpos;
	ave_abd = 0;
	dev_abd = 0;

	locate_ascending_position();
	locate_descending_position();

	check_empty();
	if(empty == false) estimate_abundance();
}

region::region(const region &r)
	:lpos(r.lpos), rpos(r.rpos), imap(r.imap), ltype(r.ltype), rtype(r.rtype)
{
	asc_pos = r.asc_pos;
	desc_pos = r.desc_pos;
	ave_abd = r.ave_abd;
	dev_abd = r.dev_abd;
	empty = r.empty;
}

region& region::operator=(const region &r)
{
	lpos = r.lpos;
	rpos = r.rpos;
	imap = r.imap;
	ltype = r.ltype;
	rtype = r.rtype;
	asc_pos = r.asc_pos;
	desc_pos = r.desc_pos;
	return *this;
}

region::~region()
{}

int region::print(int index)
{
	char em = empty ? 'T' : 'F';
	char cl = lpos < asc_pos ? 'T' : 'F';
	char cr = rpos > desc_pos ? 'T' : 'F';
	printf("region %d: [%d-%d), type = (%d, %d), empty = %c, check = (%c, %c), core = [%d-%d), len = %d, core-len = %d, ave-abd = %.1lf, std-abd = %.1lf\n",
			index, lpos, rpos, ltype, rtype, em, cl, cr, asc_pos, desc_pos, rpos - lpos, desc_pos - asc_pos, ave_abd, dev_abd);
	return 0;
}

int region::check_empty()
{
	empty = false;
	if(ltype == RIGHT_SPLICE) return 0;
	if(rtype == LEFT_SPLICE) return 0;

	/*
	if(desc_pos <= asc_pos)
	{
		empty = true;
		return 0;
	}
	*/

	int n = (rpos - lpos < num_sample_positions) ? (rpos - lpos) : num_sample_positions;
	int t = (rpos - lpos) / n;
	int r = (rpos - lpos) / t;

	int x = compute_coverage(*imap, lpos, rpos, t);
	if(1.0 * x / r >= min_region_coverage) return 0;

	int s = maximum_overlap(*imap, lpos, rpos, t);
	if(s >= min_max_region_overlap) return 0;

	empty = true;
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

int region::estimate_abundance()
{
	/*
	if(desc_pos <= asc_pos)
	{
		ave_abd = 0;
		dev_abd = 0;
		return 0;
	}
	*/

	int n = (rpos - lpos < num_sample_positions) ? (rpos - lpos) : num_sample_positions;
	int t = (rpos - lpos) / n;

	vector<int> v;
	int sum = 0;
	for(int32_t p = lpos; p < rpos; p += t)
	{
		int a = compute_overlap(*imap, p);
		sum += a;
		v.push_back(a);
	}

	ave_abd = 1.0 * sum / v.size();

	double var = 0;
	for(int i = 0; i < v.size(); i++)
	{
		var += (v[i] - ave_abd) * (v[i] - ave_abd);	
	}

	dev_abd = sqrt(var / v.size());
	return 0;
}

bool region::left_break()
{
	if(asc_pos > lpos) return true;
	else return false;
}

bool region::right_break()
{
	if(desc_pos < rpos) return true;
	else return false;
}
