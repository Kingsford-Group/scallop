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
}

region::~region()
{}

vector<partial_exon> region::build()
{
	init();
	if(empty == true) return pexons;

	build_bins();
	build_slopes();
	select_slopes();
	build_partial_exons();
	return pexons;
}

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

int region::evaluate_rectangle(int ll, int rr, double &ave, double &dev)
{
	ave = 0;
	dev = 1.0;

	if(empty == true) return 0;

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

		//printf("pos = [%d, %d), value = %d\n", lower(it->first), upper(it->first), it->second);

		if(it == rit) break;
	}

	dev = sqrt(var / (rr - ll));
	if(dev < 1.0) dev = 1.0;

	return 0;
}

int region::evaluate_triangle(int ll, int rr, double &ave, double &dev)
{
	ave = 0;
	dev = 1.0;

	if(empty == true) return 0;

	SIMI lit, rit;
	tie(lit, rit) = locate_boundary_iterators(*imap, ll, rr);

	if(lit == imap->end()) return 0;
	if(rit == imap->end()) return 0;

	vector<double> xv;
	vector<double> yv;
	double xm = 0;
	double ym = 0;
	for(SIMI it = lit; ; it++)
	{
		assert(upper(it->first) > lower(it->first));
		double xi = (lower(it->first) + upper(it->first)) / 2.0;
		double yi = it->second;
		xv.push_back(xi);
		yv.push_back(yi);
		xm += xi;
		ym += yi;
		if(it == rit) break;
	}

	xm /= xv.size();
	ym /= yv.size();

	double f1 = 0;
	double f2 = 0;
	for(int i = 0; i < xv.size(); i++)
	{
		f1 += (xv[i] - xm) * (yv[i] - ym);
		f2 += (xv[i] - xm) * (xv[i] - xm);
	}

	double b1 = f1 / f2;
	double b0 = ym - b1 * xm;

	double a1 = b1 * rr + b0;
	double a0 = b1 * ll + b0;
	ave = (a1 > a0) ? a1 : a0;

	double var = 0;
	for(SIMI it = lit; ; it++)
	{
		assert(upper(it->first) > lower(it->first));
		double xi = (upper(it->first) + lower(it->first)) / 2.0;
		double yi = b1 * xi + b0;
		var += (it->second - yi) * (it->second - yi) * (upper(it->first) - lower(it->first));
		if(it == rit) break;
	}

	dev = sqrt(var / (rr - ll));
	if(dev < 1.0) dev = 1.0;

	return 0;
}

double region::compute_deviation(const split_interval_map &sim, const vector<slope> &ss)
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
		evaluate_rectangle(l, r, a, d);
		var += d * d * (r - l);
	}

	for(int i = 0; i < ss.size(); i++)
	{
		const slope &s = ss[i];
		var += s.dev * s.dev * (s.rpos - s.lpos);
		len += (s.rpos - s.lpos);
	}
	return sqrt(var / len);
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
	
	for(int i = 0; i < bins.size(); i++)
	{
		bins[i] = bins[i] / bsize;
	}

	return 0;
}

int region::build_slopes()
{
	for(int i = slope_min_bin_num; i <= slope_max_bin_num; i += 3)
	{
		build_slopes(i);
	}
	return 0;
}

int region::build_slopes(int bin_num)
{
	if(bins.size() < bin_num) return 0;

	assert(bin_num % 3 == 0);
	int d = bin_num / 3;

	for(int i = 0; i < bins.size() - bin_num; i++)
	{
		int xo = 0, yo = 0, zo = 0;
		for(int k = 0; k < d; k++)
		{
			xo += bins[i + k + 0 * d];
			yo += bins[i + k + 1 * d];
			zo += bins[i + k + 2 * d];
		}
		int xx = xo * slope_bin_size / average_read_length;
		int yy = yo * slope_bin_size / average_read_length;
		int zz = zo * slope_bin_size / average_read_length;

		int sxy = compute_binomial_score(xx + yy, 0.5, yy);
		int syz = compute_binomial_score(yy + zz, 0.5, zz);

		int score5 = (sxy < syz) ? sxy : syz;

		if(score5 > slope_min_score)
		{
			slope s5(SLOPE5END, i, i + bin_num, score5);
			extend_slope(s5);
			seeds.push_back(s5);
		}

		int syx = compute_binomial_score(xx + yy, 0.5, xx);
		int szy = compute_binomial_score(yy + zz, 0.5, yy);

		int score3 = (syx < szy) ? syx : szy;

		if(score3 > slope_min_score)
		{
			slope s3(SLOPE3END, i, i + bin_num, score3);
			extend_slope(s3);
			seeds.push_back(s3);
		}
	}

	return 0;
}

int region::extend_slope(slope &s)
{
	for(int i = s.lbin - 1; i >= 0; i--)
	{
		if(s.type == SLOPE5END && bins[i] < bins[i + 1]) continue;
		if(s.type == SLOPE3END && bins[i] > bins[i + 1]) continue;
		s.lbin = i + 1;
		break;
	}
	for(int i = s.rbin; i < bins.size(); i++)
	{
		if(s.type == SLOPE5END && bins[i] > bins[i - 1]) continue;
		if(s.type == SLOPE3END && bins[i] < bins[i - 1]) continue;
		s.rbin = i;
		break;
	}

	s.lpos = lcore + s.lbin * slope_bin_size;

	if(s.rbin == bins.size()) s.rpos = rcore;
	else s.rpos = lcore + s.rbin * slope_bin_size;

	if(s.lpos <= lcore + slope_flexible_size) s.lpos = lcore;
	if(s.rpos >= rcore - slope_flexible_size) s.rpos = rcore;

	if(s.lpos == lcore || s.rpos == rcore) s.flag = SLOPE_MARGIN;
	else s.flag = SLOPE_MIDDLE;

	double ave, dev;
	evaluate_triangle(s.lpos, s.rpos, s.ave, s.dev);

	return 0;
}

int region::select_slopes()
{
	split_interval_map sim;
	sim += make_pair(ROI(lcore, rcore), 1);

	double dev = compute_deviation(sim, slopes);

	std::sort(seeds.begin(), seeds.end(), compare_slope_score);
	for(int i = 0; i < seeds.size(); i++)
	{
		slope &x = seeds[i];

		bool b = true;
		for(int k = 0; k < slopes.size(); k++)
		{
			slope &y = slopes[k];
			int d = x.distance(y);
			if(d <= slope_min_distance) b = false;
			if(b == false) break;
		}
		if(b == false) continue;

		x.print(i);

		slopes.push_back(x);
		sim -= make_pair(ROI(x.lpos, x.rpos), 1);

		double d = compute_deviation(sim, slopes);

		if(d <= dev - dev * slope_acceptance_dev_decrease)
		{
			dev = d;
		}
		else
		{
			slopes.pop_back();
			break;
		}
	}

	return 0;
}

int region::build_partial_exons()
{
	std::sort(slopes.begin(), slopes.end(), compare_slope_pos);

	int32_t lexon = lcore;
	int32_t rexon = -1;
	int lltype = ltype;
	int rrtype = -1;
	for(int i = 0; i < slopes.size(); i++)
	{
		slope &s = slopes[i];

		//s.print(i);

		if(s.type == SLOPE5END)
		{
			rexon = s.lpos;
			rrtype = START_BOUNDARY;

			if(lexon < rexon)
			{
				partial_exon pe(lexon, rexon, lltype, rrtype);
				evaluate_rectangle(pe.lpos, pe.rpos, pe.ave_abd, pe.dev_abd);
				pexons.push_back(pe);
			}

			lexon = rexon;
			lltype = rrtype;
			rexon = s.rpos;
			rrtype = (rexon == rcore) ? rtype : MIDDLE_CUT;

			partial_exon pe(lexon, rexon, lltype, rrtype);
			evaluate_rectangle(pe.lpos, pe.rpos, pe.ave_abd, pe.dev_abd);
			pexons.push_back(pe);

			lexon = rexon;
			lltype = rrtype;
		}
		else if(s.type == SLOPE3END)
		{
			rexon = s.lpos;
			rrtype = (rexon == lcore) ? ltype : MIDDLE_CUT;

			if(lexon < rexon)
			{
				partial_exon pe(lexon, rexon, lltype, rrtype);
				evaluate_rectangle(pe.lpos, pe.rpos, pe.ave_abd, pe.dev_abd);
				pexons.push_back(pe);
			}

			lexon = rexon;
			lltype = rrtype;
			rexon = s.rpos;
			rrtype = (rexon == rcore) ? rtype : END_BOUNDARY;

			partial_exon pe(lexon, rexon, lltype, rrtype);
			evaluate_rectangle(pe.lpos, pe.rpos, pe.ave_abd, pe.dev_abd);
			pexons.push_back(pe);

			lexon = rexon;
			lltype = rrtype;
		}
	}

	rexon = rcore;
	rrtype = rtype;

	if(lexon < rexon)
	{
		partial_exon pe(lexon, rexon, lltype, rrtype);
		evaluate_rectangle(pe.lpos, pe.rpos, pe.ave_abd, pe.dev_abd);
		pexons.push_back(pe);
	}

	return 0;
}

int region::print(int index) const
{
	int32_t lc = compute_overlap(*imap, lcore);
	int32_t rc = compute_overlap(*imap, rcore - 1);
	printf("region %d: empty = %c, pos = [%d, %d), core = [%d, %d), bins = %lu, coverage = (%d, %d)\n", 
			index, empty ? 'T' : 'F', lpos, rpos, lcore, rcore, bins.size(), lc, rc);


	/*
	for(int i = 0; i < seeds.size(); i++)
	{
		seeds[i].print(i);
	}
	*/
	for(int i = 0; i < slopes.size(); i++)
	{
		slopes[i].print(i);
	}

	return 0;
	for(int i = 0; i < pexons.size(); i++)
	{
		pexons[i].print(i);
	}
	printf("\n");
	return 0;
}


