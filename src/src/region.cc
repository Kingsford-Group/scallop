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
	if(bins.size() < slope_min_bin_num) return 0;

	// left part
	int mbin = slope_min_bin_num / 3;
	int nbin = slope_std_bin_num / 3;

	int lbin = 0 * mbin;
	int rbin = 3 * (mbin - 1);
	int xo = 0;
	int yo = 0;
	int zo = 0;
	int d = 0;
	for(d = 0; d < mbin - 1; d++)
	{
		xo += bins[0 * (mbin - 1) + d];
		yo += bins[1 * (mbin - 1) + d];
		zo += bins[2 * (mbin - 1) + d];
	}
	assert(d == mbin - 1);

	seeds.clear();
	for(d = mbin - 1; d < nbin && rbin + 3 <= bins.size(); d++)
	{
		xo += bins[lbin + 1 * d + 0];
		yo -= bins[lbin + 1 * d + 0];
		yo += bins[lbin + 2 * d + 0];
		yo += bins[lbin + 2 * d + 1];
		zo -= bins[lbin + 2 * d + 0];
		zo -= bins[lbin + 2 * d + 1];
		zo += bins[lbin + 3 * d + 0];
		zo += bins[lbin + 3 * d + 1];
		zo += bins[lbin + 3 * d + 2];

		rbin += 3;

		int xx = xo / average_read_length;
		int yy = yo / average_read_length;
		int zz = zo / average_read_length;

		int sxy = compute_binomial_score(xx + yy, 0.5, yy);
		int syz = compute_binomial_score(yy + zz, 0.5, zz);

		slope s5(SLOPE5END, lbin, rbin, (sxy < syz) ? sxy : syz);
		if(s5.score > min_slope_score) seeds.push_back(s5);

		int syx = compute_binomial_score(xx + yy, 0.5, xx);
		int szy = compute_binomial_score(yy + zz, 0.5, yy);
		slope s3(SLOPE3END, lbin, rbin, (syx < szy) ? syx : szy);
		if(s3.score > min_slope_score) seeds.push_back(s3);
	}

	// middle part
	while(rbin + 1 <= bins.size())
	{
		xo -= bins[lbin + 0 * d];
		xo += bins[lbin + 1 * d];
		yo -= bins[lbin + 1 * d];
		yo += bins[lbin + 2 * d];
		zo -= bins[lbin + 2 * d];
		zo += bins[lbin + 3 * d];

		lbin++;
		rbin++;

		int xx = xo / average_read_length;
		int yy = yo / average_read_length;
		int zz = zo / average_read_length;

		int sxy = compute_binomial_score(xx + yy, 0.5, yy);
		int syz = compute_binomial_score(yy + zz, 0.5, zz);

		slope s5(SLOPE5END, lbin, rbin, (sxy < syz) ? sxy : syz);
		if(s5.score > min_slope_score) seeds.push_back(s5);

		int syx = compute_binomial_score(xx + yy, 0.5, xx);
		int szy = compute_binomial_score(yy + zz, 0.5, yy);
		slope s3(SLOPE3END, lbin, rbin, (syx < szy) ? syx : szy);
		if(s3.score > min_slope_score) seeds.push_back(s3);
	}

	assert(rbin == bins.size());

	// right part
	for(; d > mbin; d--) 
	{
		xo -= bins[lbin + 0 * d + 0];
		xo -= bins[lbin + 0 * d + 1];
		xo -= bins[lbin + 0 * d + 2];
		xo += bins[lbin + 1 * d + 0];
		xo += bins[lbin + 1 * d + 1];

		yo -= bins[lbin + 1 * d + 0];
		yo -= bins[lbin + 1 * d + 1];
		yo += bins[lbin + 2 * d + 0];

		zo -= bins[lbin + 2 * d + 0];

		lbin += 3;

		int xx = xo / average_read_length;
		int yy = yo / average_read_length;
		int zz = zo / average_read_length;

		int sxy = compute_binomial_score(xx + yy, 0.5, yy);
		int syz = compute_binomial_score(yy + zz, 0.5, zz);

		slope s5(SLOPE5END, lbin, rbin, (sxy < syz) ? sxy : syz);
		if(s5.score > min_slope_score) seeds.push_back(s5);

		int syx = compute_binomial_score(xx + yy, 0.5, xx);
		int szy = compute_binomial_score(yy + zz, 0.5, yy);
		slope s3(SLOPE3END, lbin, rbin, (syx < szy) ? syx : szy);
		if(s3.score > min_slope_score) seeds.push_back(s3);
	}

	for(int i = 0; i < seeds.size(); i++)
	{
		slope &s = seeds[i];
		s.lpos = lcore + s.lbin * slope_bin_size;

		if(s.rbin == bins.size()) s.rpos = rcore;
		else s.rpos = lcore + s.rbin * slope_bin_size;

		if(s.lpos == lcore || s.rpos == rcore) s.flag = SLOPE_MARGIN;
		else s.flag = SLOPE_MIDDLE;

		double ave, dev;
		evaluate_triangle(s.lpos, s.rpos, s.ave, s.dev);
	}

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
		//evaluate_triangle(x.lpos, x.rpos, x.ave, x.dev);	

		bool b = true;
		for(int k = 0; k < slopes.size(); k++)
		{
			slope &y = slopes[k];
			int d = x.distance(y);
			if(d <= slope_min_distance) b = false;
			if(b == false) break;
		}
		if(b == false) continue;

		slopes.push_back(x);
		sim -= make_pair(ROI(x.lpos, x.rpos), 1);

		double d = compute_deviation(sim, slopes);

		if(d <= dev * 0.9)
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
	int32_t ppos = lcore;
	double pvar = 0.0;
	double psum = 0.0;
	int lltype = ltype;
	for(int i = 0; i < slopes.size(); i++)
	{
		slope &s = slopes[i];
		if(s.type == SLOPE5END)
		{
			int rexon = s.lpos;
			int rrtype = START_BOUNDARY;

			if(rexon > ppos)
			{
				double ave, dev;
				evaluate_rectangle(ppos, rexon, ave, dev);
				psum += ave * (rexon - ppos);
				pvar += dev * dev * (rexon - ppos);
			}

			if(lexon < rexon)
			{
				double ave = psum / (rexon - lexon);
				double dev = sqrt(pvar / (rexon - lexon));
				partial_exon pe(lexon, rexon, lltype, rrtype);
				pe.ave_abd = ave;
				pe.dev_abd = dev;
				pexons.push_back(pe);
			}

			lexon = s.lpos;
			ppos = s.rpos;
			lltype = START_BOUNDARY;
			psum = s.ave * (s.rpos - s.lpos);
			pvar = s.dev * s.dev * (s.rpos - s.lpos);
		}
		else if(s.type == SLOPE3END)
		{
			int rexon = s.rpos;
			int rrtype = END_BOUNDARY;

			if(s.lpos > ppos)
			{
				double ave, dev;
				evaluate_rectangle(ppos, s.lpos, ave, dev);
				psum += ave * (s.lpos - ppos);
				pvar += dev * dev * (s.lpos - ppos);
			}

			psum += s.ave * (s.rpos - s.lpos);
			pvar += s.dev * s.dev * (s.rpos - s.lpos);

			assert(lexon < rexon);

			double ave = psum / (rexon - lexon);
			double dev = sqrt(pvar / (rexon - lexon));
			partial_exon pe(lexon, rexon, lltype, rrtype);
			pe.ave_abd = ave;
			pe.dev_abd = dev;
			pexons.push_back(pe);

			lexon = s.rpos;
			ppos = s.rpos;
			lltype = END_BOUNDARY;
			psum = 0;
			pvar = 0;
		}
	}

	int rexon = rcore;
	int rrtype = rtype;

	if(rexon > ppos)
	{
		double ave, dev;
		//printf("-------------------\n");
		//printf("rectangle [%d, %d), ave = %.2lf, dev = %.2lf\n", ppos, rexon, ave, dev);
		//printf("===================\n");
		evaluate_rectangle(ppos, rexon, ave, dev);
		psum += ave * (rexon - ppos);
		pvar += dev * dev * (rexon - ppos);
	}

	if(lexon < rexon)
	{
		double ave = psum / (rexon - lexon);
		double dev = sqrt(pvar / (rexon - lexon));
		partial_exon pe(lexon, rexon, lltype, rrtype);
		pe.ave_abd = ave;
		pe.dev_abd = dev;
		pexons.push_back(pe);
	}

	return 0;
}

int region::print(int index) const
{
	printf("region %d: empty = %c, pos = [%d, %d), core = [%d, %d), bins = %lu\n", 
			index, empty ? 'T' : 'F', lpos, rpos, lcore, rcore, bins.size());
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
	for(int i = 0; i < pexons.size(); i++)
	{
		pexons[i].print(i);
	}
	printf("\n");
	return 0;
}


