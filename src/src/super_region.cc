#include "super_region.h"
#include "config.h"
#include "util.h"
#include "binomial.h"
#include <algorithm>

using namespace std;

super_region::super_region(const split_interval_map *_imap)
	:imap(_imap)
{
}

super_region::~super_region()
{}

int super_region::clear()
{
	bins.clear();
	seeds.clear();
	slopes.clear();
	pexons.clear();
	return 0;
}

int super_region::build()
{
	build_bins();
	build_slopes();
	select_slopes(0, bins.size(), 0);
	adjust_coverage();
	assign_boundaries();
	build_partial_exons();
	return 0;
}

int super_region::add_region(const region &r)
{
	regions.push_back(r);
	return 0;
}

int super_region::build_bins()
{
	bins.clear();
	for(int i = 0; i < regions.size(); i++)
	{
		bins.insert(bins.end(), regions[i].bins.begin(), regions[i].bins.end());
	}
	return 0;
}

int super_region::build_slopes()
{
	int bin_num = slope_min_bin_num;

	if(bins.size() < bin_num) return 0;

	assert(bin_num % 4 == 0);
	int d = bin_num / 4;

	for(int i = 0; i < bins.size() - bin_num; i++)
	{
		int wo = 0, xo = 0, yo = 0, zo = 0;
		for(int k = 0; k < d; k++)
		{
			wo += bins[i + k + 0 * d];
			xo += bins[i + k + 1 * d];
			yo += bins[i + k + 2 * d];
			zo += bins[i + k + 3 * d];
		}
		int ww = wo * slope_bin_size / average_read_length;
		int xx = xo * slope_bin_size / average_read_length;
		int yy = yo * slope_bin_size / average_read_length;
		int zz = zo * slope_bin_size / average_read_length;

		int swx = compute_binomial_score(ww + xx, 0.5, xx);
		int sxy = compute_binomial_score(xx + yy, 0.5, yy);
		int syz = compute_binomial_score(yy + zz, 0.5, zz);

		int score5 = (swx < sxy && swx < syz) ? swx : (sxy < syz ? sxy : syz);

		if(score5 > slope_min_score)
		{
			slope s5(SLOPE5END, i, i + bin_num, score5);
			extend_slope(s5);
			evaluate_slope(s5);
			seeds.push_back(s5);
		}

		int sxw = compute_binomial_score(ww + xx, 0.5, ww);
		int syx = compute_binomial_score(xx + yy, 0.5, xx);
		int szy = compute_binomial_score(yy + zz, 0.5, yy);

		int score3 = (sxw < syx && sxw < szy) ? sxw : (syx < szy ? syx : szy);

		if(score3 > slope_min_score)
		{
			slope s3(SLOPE3END, i, i + bin_num, score3);
			extend_slope(s3);
			evaluate_slope(s3);
			seeds.push_back(s3);
		}
	}

	std::sort(seeds.begin(), seeds.end(), compare_slope_score);

	//for(int i = 0; i < seeds.size(); i++) seeds[i].print(i);

	return 0;
}

int super_region::extend_slope(slope &s)
{
	int m = s.rbin - s.lbin;
	assert(m % 4 == 0);
	m = m / 4;

	int tail_coverage = 8;

	// extend to the left 
	int p = 0;
	for(int i = s.lbin; i < s.lbin + m; i++) p += bins[i];
	for(int i = s.lbin - m; i >= 0; i -= m)
	{
		int q = 0;
		for(int j = i; j < i + m; j++) q += bins[j];
		int score = 0;
		if(s.type == SLOPE5END) score = compute_binomial_score(p + q, 0.5, p);
		if(s.type == SLOPE3END) score = compute_binomial_score(p + q, 0.5, q);
		if(score >= slope_extend_score) s.lbin -= m;
		else break;
		p = q;
	}

	// extend to the right
	p = 0;
	for(int i = s.rbin - m; i < s.rbin; i++) p += bins[i];
	for(int i = s.rbin; i + m <= bins.size(); i += m)
	{
		int q = 0;
		for(int j = i; j < i + m; j++) q += bins[j];
		int score = 0;
		if(s.type == SLOPE5END) score = compute_binomial_score(p + q, 0.5, q);
		if(s.type == SLOPE3END) score = compute_binomial_score(p + q, 0.5, p);
		if(score >= slope_extend_score) s.rbin += m;
		else break;
		p = q;
	}

	// further small extend
	for(int i = s.lbin - 1; i >= 0; i--)
	{
		if(s.type == SLOPE5END && bins[i] < bins[i + 1]) s.lbin = i;
		else if(s.type == SLOPE3END && bins[i] > bins[i + 1]) s.lbin = i;
		else if(s.type == SLOPE5END && bins[i] <= tail_coverage) s.lbin = i;
		else break;
	}
	for(int i = s.rbin; i < bins.size(); i++)
	{
		if(s.type == SLOPE5END && bins[i] > bins[i - 1]) s.rbin = i + 1;
		else if(s.type == SLOPE3END && bins[i] < bins[i - 1]) s.rbin = i + 1;
		else if(s.type == SLOPE3END && bins[i] <= tail_coverage) s.rbin = i + 1;
		else break;
	}

	// adjust the starting/end positions
	int f = slope_flexible_bin_num;
	int xi, xb;
	int yi, yb;
	locate_bin(s.lbin, xi, xb);
	locate_bin(s.rbin - 1, yi, yb);
	int nx = regions[xi].bins.size();
	int ny = regions[yi].bins.size();
	if(nx >= 2 * f && xb <= f) s.lbin -= xb;
	if(nx >= 2 * f && xb + f >= nx) s.lbin += (nx - xb);
	if(ny >= 2 * f && yb < f) s.rbin -= (yb + 1);
	if(ny >= 2 * f && yb + f >= ny - 1) s.rbin += (ny - yb - 1);

	// handle start and end boundary
	locate_bin(s.lbin, xi, xb);
	locate_bin(s.rbin - 1, yi, yb);
	nx = regions[xi].bins.size();
	ny = regions[yi].bins.size();
	//&& bins[s.lbin] <= 2 * tail_coverage
	if(s.type == SLOPE5END && regions[xi].ltype == START_BOUNDARY  && xb < average_slope_length / slope_bin_size)
	{
		s.lbin -= xb;
	}
	if(s.type == SLOPE3END && regions[yi].rtype == END_BOUNDARY  && ny - yb - 1 < average_slope_length / slope_bin_size)
	{
		s.rbin += (ny - yb - 1);
	}

	assert(s.lbin < s.rbin);
	assert(s.lbin >= 0 && s.lbin < bins.size());
	assert(s.rbin > 0 && s.rbin <= bins.size());

	locate_bin(s.lbin, xi, xb);
	locate_bin(s.rbin - 1, yi, yb);

	s.lpos = regions[xi].get_location(xb);
	s.rpos = regions[yi].get_location(yb + 1);

	return 0;
}

int super_region::evaluate_slope(slope &s)
{
	compute_mean_dev(bins, s.lbin, s.rbin, s.ave, s.dev);
	return 0;
}

int super_region::locate_bin(int x, int &xi, int &xb)
{
	xi = xb = -1;
	int sum = 0;
	for(int i = 0; i < regions.size(); i++)
	{
		region &r = regions[i];
		if(x >= sum && x < sum + r.bins.size())
		{
			xi = i;
			xb = x - sum;
			break;
		}
		sum += r.bins.size();
	}

	//printf("locate %d to (%d, %d), total = %lu\n", x, xi, xb, bins.size());
	assert(xi >= 0);
	assert(xb >= 0 && xb < regions[xi].bins.size());

	return 0;
}

int super_region::select_slopes(int si, int ti, int ss)
{
	double ave, dev;
	compute_mean_dev(bins, si, ti, ave, dev);

	int k = -1;
	for(int i = ss; i < seeds.size(); i++)
	{
		slope &x = seeds[i];
		if(x.lbin < si) continue;
		if(x.rbin > ti) continue;
		k = i;
		break;
	}

	if(k == -1) return 0;

	slope &x = seeds[k];

	double ave1, dev1;
	double ave2, dev2;

	int bm = average_slope_length / slope_bin_size;
	int ssi = (si > x.lbin - bm) ? si : x.lbin - bm;
	int tti = (ti < x.rbin + bm) ? ti : x.rbin + bm;

	compute_mean_dev(bins, ssi, x.lbin, ave1, dev1);
	compute_mean_dev(bins, x.rbin, tti, ave2, dev2);

	bool b1 = (x.lbin - ssi >= 5) && (x.type == SLOPE5END) && (x.ave <= ave1 + dev1 * slope_acceptance_sigma);
	bool b2 = (x.lbin - ssi >= 5) && (x.type == SLOPE3END) && (x.ave >= ave1 - dev1 * slope_acceptance_sigma);
	bool b3 = (tti - x.rbin >= 5) && (x.type == SLOPE5END) && (x.ave <= ave2 - dev2 * slope_acceptance_sigma);
	bool b4 = (tti - x.rbin >= 5) && (x.type == SLOPE3END) && (x.ave >= ave2 + dev2 * slope_acceptance_sigma);

	//x.print(k);
	//printf("select flag = (%c, %c): (%d, %d, %.2lf, %.2lf) -> (%d, %d, %.2lf, %.2lf) + (%d, %d, %.2lf, %.2lf)\n", 
	//	b1 ? 'F' : 'T', b2 ? 'F' : 'T', si, ti, ave, dev, si, x.lbin, ave1, dev1, x.rbin, ti, ave2, dev2);

	if(b1 == false && b2 == false && b3 == false && b4 == false)
	{
		slopes.push_back(x);
		select_slopes(si, x.lbin, 0);
		select_slopes(x.rbin, ti, 0);
	}
	else
	{
		select_slopes(si, ti, k + 1);
	}

	return 0;
}

int super_region::adjust_coverage()
{
	for(int i = 0; i < slopes.size(); i++)
	{
		slope &s = slopes[i];

		adjust_coverage(s.lbin, s.rbin);
		evaluate_slope(s);

		s.print(i);

		int xi, xb;
		int yi, yb;
		locate_bin(s.lbin, xi, xb);
		locate_bin(s.rbin - 1, yi, yb);

		int k = s.lbin;
		for(int ii = xi; ii <= yi; ii++)
		{
			int xx = (ii == xi) ? xb : 0;
			int yy = (ii == yi) ? yb + 1 : regions[ii].bins.size();
			for(int j = xx; j < yy; j++)
			{
				regions[ii].bins[j] = bins[k++];
				regions[ii].adjust[j] = true;
			}
		}
	}
	return 0;
}

int super_region::adjust_coverage(int si, int ti)
{
	if(ti <= si) return 0;

	int n = (ti - si) / 2;
	int s1 = 0;
	int s2 = 0;
	for(int i = 0; i < n; i++)
	{
		s1 += bins[si + i * 2 + 0];
		s2 += bins[si + i * 2 + 1];
	}
	double d = (s2 - s1) * 1.0 / n;

	for(int i = si; i < ti; i++)
	{
		if(d >= 0) bins[i] += d * (ti - i - 1);
		else bins[i] += (0.0 - d) * (i - si);
	}
	return 0;
}

int super_region::build_partial_exons()
{
	pexons.clear();
	for(int i = 0; i < regions.size(); i++)
	{
		vector<partial_exon> v;
		regions[i].build_partial_exons(v);
		pexons.insert(pexons.end(), v.begin(), v.end());
	}
	return 0;
}

int super_region::assign_boundaries()
{
	for(int i = 0; i < slopes.size(); i++)
	{
		slope &s = slopes[i];
		int xi, xb;
		int yi, yb;
		if(s.type == SLOPE5END)
		{
			locate_bin(s.lbin, xi, xb);
			locate_bin(s.rbin - 1, yi, yb);
			regions[xi].add_boundary(xb, START_BOUNDARY);
			regions[yi].add_boundary(yb + 1, MIDDLE_CUT);
		}
		else if(s.type == SLOPE3END)
		{
			locate_bin(s.lbin, xi, xb);
			locate_bin(s.rbin - 1, yi, yb);
			regions[xi].add_boundary(xb, MIDDLE_CUT);
			regions[yi].add_boundary(yb + 1, END_BOUNDARY);
		}
	}
	return 0;
}

int super_region::size() const
{
	return regions.size();
}

int super_region::print(int index) const
{
	printf("super region %d: \n", index);
	for(int i = 0; i < regions.size(); i++)
	{
		printf(" ");
		regions[i].print(i);
	}

	/*
	for(int i = 0; i < bins.size(); i++)
	{
		printf("super region %d bin %d = %d\n", index, i, bins[i]);
	}
	printf("\n");

	for(int i = 0; i < slopes.size(); i++)
	{
		printf(" ");
		slopes[i].print(i);
	}
	*/
	return 0;
}
