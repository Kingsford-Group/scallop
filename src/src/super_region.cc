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
	assign_coverage();
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
	for(int i = slope_min_bin_num; i <= slope_max_bin_num; i += 4)
	{
		build_slopes(i);
	}
	std::sort(seeds.begin(), seeds.end(), compare_slope_score);
	
	//for(int i = 0; i < seeds.size(); i++) seeds[i].print(i);

	return 0;
}

int super_region::build_slopes(int bin_num)
{
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

	return 0;
}

int super_region::extend_slope(slope &s)
{
	// extend the given slope
	for(int i = s.lbin - 1; i >= 0; i--)
	{
		if(s.type == SLOPE5END && bins[i] < bins[i + 1]) s.lbin = i;
		if(s.type == SLOPE3END && bins[i] > bins[i + 1]) s.lbin = i;
		break;
	}
	for(int i = s.rbin; i < bins.size(); i++)
	{
		if(s.type == SLOPE5END && bins[i] > bins[i - 1]) s.rbin = i + 1;
		if(s.type == SLOPE3END && bins[i] < bins[i - 1]) s.rbin = i + 1;
		break;
	}

	int f = slope_flexible_bin_num;
	int xi, xb;

	locate_bin(s.lbin, xi, xb);
	int n = regions[xi].bins.size();
	if(n >= 2 * f && xb <= f) s.lbin -= xb;
	if(n >= 2 * f && xb + f >= n) s.lbin += (n - xb);

	locate_bin(s.rbin - 1, xi, xb);
	n = regions[xi].bins.size();
	if(n >= 2 * f && xb < f) s.rbin -= (xb + 1);
	if(n >= 2 * f && xb + f >= n - 1) s.rbin += (n - xb - 1);

	assert(s.lbin < s.rbin);
	assert(s.lbin >= 0 && s.lbin < bins.size());
	assert(s.rbin > 0 && s.rbin <= bins.size());

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

	int ssi = (si > x.lbin - slope_max_bin_num) ? si : x.lbin - slope_max_bin_num;
	int tti = (ti < x.rbin + slope_max_bin_num) ? ti : x.rbin + slope_max_bin_num;

	compute_mean_dev(bins, ssi, x.lbin, ave1, dev1);
	compute_mean_dev(bins, x.rbin, tti, ave2, dev2);

	bool b1 = (x.lbin - ssi >= 5) && (fabs(ave1 - x.ave) < dev1 * slope_acceptance_sigma);
	bool b2 = (tti - x.rbin >= 5) && (fabs(ave2 - x.ave) < dev2 * slope_acceptance_sigma);

	//x.print(k);
	//printf("select (%d, %d, %.2lf, %.2lf) -> (%d, %d, %.2lf, %.2lf) + (%d, %d, %.2lf, %.2lf)\n", si, ti, ave, dev, si, x.lbin, ave1, dev1, x.rbin, ti, ave2, dev2);

	if(b1 == false && b2 == false)
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
		s.print(i);

		int xi, xb;
		if(s.type == SLOPE5END) locate_bin(s.lbin, xi, xb);
		if(s.type == SLOPE3END) locate_bin(s.rbin - 1, xi, xb);
		if(s.type == SLOPE5END) regions[xi].add_boundary(xb, s.type);
		if(s.type == SLOPE3END) regions[xi].add_boundary(xb + 1, s.type);
	}
	return 0;
}

int super_region::assign_coverage()
{
	int k = 0;
	for(int i = 0; i < regions.size(); i++)
	{
		region &r = regions[i];
		for(int j = 0; j < r.bins.size(); j++)
		{
			r.bins[j] = bins[k];
			k++;
		}
	}
	assert(k == bins.size());
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
