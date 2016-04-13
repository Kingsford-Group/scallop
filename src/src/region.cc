#include "region.h"
#include "config.h"
#include "util.h"
#include "binomial.h"
#include <algorithm>

using namespace std;

region::region(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype, const split_interval_map *_imap)
	:lpos(_lpos), rpos(_rpos), imap(_imap), ltype(_ltype), rtype(_rtype)
{
	ave_abd = 0;
	dev_abd = 1;
	lcore = lpos;
	rcore = rpos;
}

region::~region()
{}

int region::build()
{
	init_core_empty();
	estimate_abundance();
	compute_bin_abundances();
	compute_slopes();
	select_slopes();
	return 0;
}

int region::init_core_empty()
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
	return estimate_abundance(lcore, rcore, ave_abd, dev_abd);
}

double region::compute_deviation(const split_interval_map &sim)
{
	double var = 0;
	int len = 0;
	for(SIMI it = sim.begin(); it != sim.end(); it++)
	{
		int l = lower(it->first);
		int r = upper(it->first);
		len += (r - l);
		int w = it->second;
		if(w <= 0) continue;

		double a, d;
		estimate_abundance(l, r, a, d);
		var += d * d * (r - l);
	}
	return sqrt(var / len);
}

int region::estimate_abundance(int ll, int rr, double &ave, double &dev)
{
	if(empty == true) return 0;

	ave = 0;
	dev = 0;

	SIMI lit, rit;
	tie(lit, rit) = locate_boundary_iterators(*imap, ll, rr);

	if(lit == imap->end()) return 0;
	if(rit == imap->end()) return 0;

	ave = 1.0 * compute_sum_overlap(*imap, lit, rit) / (rr - ll);

	double var = 0;
	for(SIMI it = lit; ; it++)
	{
		assert(upper(it->first) > lower(it->first));
		var += (it->second - ave) * (it->second - ave) * (upper(it->first) - lower(it->first));
		if(it == rit) break;
	}

	dev = sqrt(var / (rr - ll));

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

string region::label() const
{
	string l = tostring(lpos % 100000);
	string r = tostring(rpos % 100000);
	return (l + "-" + r);
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

int region::print_boundaries(int index) const
{
	int bsize = region_bin_size;
	printf("region %d slopes:\n", index);
	for(int i = 0; i < slopes.size(); i++)
	{
		//if(slopes[i].score <= 40) continue;
		const slope &s = slopes[i];
		printf("slope %d: type = %d, (%d - %d), score = %d, dev = %.2lf\n", i, s.type, s.lpos, s.rpos, s.score, s.dev);
	}

	/*
	printf("region %d bin abundances:\n", index);
	for(int i = 0; i < bins.size(); i++)
	{
		printf("%d + %d: %d\n", lcore + i * bsize, bsize, bins[i]);
	}
	*/
	return 0;
}

int region::compute_bin_abundances()
{
	bins.clear();
	if(empty == true) return 0;

	int32_t bsize = region_bin_size;
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
	
	for(int i = 0; i < bins.size(); i++)
	{
		bins[i] = bins[i] / bsize;
	}

	return 0;
}

int region::compute_slopes()
{
	if(bins.size() < transcript_end_bin_num) return 0;

	int nbin = transcript_end_bin_num / 3;
	int x = 0;
	int y = x + nbin;
	int z = y + nbin;
	
	int xo = 0;
	int yo = 0;
	int zo = 0;
	for(int i = 0; i < nbin; i++)
	{
		xo += bins[x + i];
		yo += bins[y + i];
		zo += bins[z + i];
	}

	int bsize = region_bin_size;
	int tsize = bsize * transcript_end_bin_num;
	slopes.clear();
	for(int i = nbin; i <= bins.size() - transcript_end_bin_num + nbin; i++)
	{
		int ll = lcore + bsize * (x + i - nbin);
		int rr = ll + tsize;
		if(i == bins.size() - transcript_end_bin_num + nbin) rr = rcore;

		int xx = xo / average_read_length;
		int yy = yo / average_read_length;
		int zz = zo / average_read_length;

		int sxy = compute_binomial_score(xx + yy, 0.5, yy);
		int syz = compute_binomial_score(yy + zz, 0.5, zz);

		slope s5(SLOPE5END, ll, rr, (sxy < syz) ? sxy : syz);
		if(s5.score > 40) slopes.push_back(s5);

		int syx = compute_binomial_score(xx + yy, 0.5, xx);
		int szy = compute_binomial_score(yy + zz, 0.5, yy);
		slope s3(SLOPE3END, ll, rr, (syx < szy) ? syx : szy);
		if(s3.score > 40) slopes.push_back(s3);

		xo -= bins[x + i - nbin];
		yo -= bins[y + i - nbin];
		zo -= bins[z + i - nbin];

		xo += bins[x + i];
		yo += bins[y + i];
		zo += bins[z + i];
	}

	std::sort(slopes.begin(), slopes.end());
	return 0;
}

int region::select_slopes()
{
	vector<slope> ss;

	split_interval_map sim;
	sim += make_pair(ROI(lcore, rcore), 1);

	double dev = compute_deviation(sim);

	for(int i = 0; i < slopes.size(); i++)
	{
		slope &x = slopes[i];
		bool b = true;
		for(int k = 0; k < ss.size(); k++)
		{
			slope &y = ss[k];
			int d = x.distance(y);
			if(d <= min_slope_distance) b = false;
			if(b == false) break;
		}
		if(b == false) continue;

		sim -= make_pair(ROI(x.lpos, x.rpos), 1);

		double d = compute_deviation(sim);
		if(d >= dev) break;

		dev = d;
		x.dev = dev;
		ss.push_back(x);
	}

	slopes = ss;
	return 0;
}
