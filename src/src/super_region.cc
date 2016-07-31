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
	select_slopes();
	build_partial_exons();
	return 0;
}

int super_region::add_region(const region &r)
{
	regions.push_back(r);
	return 0;
}

double super_region::compute_deviation(const split_interval_map &sim, const vector<slope> &ss)
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
		evaluate_rectangle(*imap, l, r, a, d);
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
	for(int i = slope_min_bin_num; i <= slope_max_bin_num; i += 3)
	{
		build_slopes(i);
	}
	return 0;
}

int super_region::build_slopes(int bin_num)
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
			seeds.push_back(s5);
		}

		int syx = compute_binomial_score(xx + yy, 0.5, xx);
		int szy = compute_binomial_score(yy + zz, 0.5, yy);

		int score3 = (syx < szy) ? syx : szy;

		if(score3 > slope_min_score)
		{
			slope s3(SLOPE3END, i, i + bin_num, score3);
			seeds.push_back(s3);
		}
	}

	return 0;
}

int super_region::select_slopes()
{
	split_interval_map sim;
	// TODO
	//sim += make_pair(ROI(lcore, rcore), 1);

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


		slopes.push_back(x);
		sim -= make_pair(ROI(x.lpos, x.rpos), 1);

		double d = compute_deviation(sim, slopes);

		x.print(i);
		//printf("previous dev = %.2lf, current dev = %.2lf\n", dev, d);

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

int super_region::build_partial_exons()
{
	return 0;
	/*
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
	*/
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

	for(int i = 0; i < slopes.size(); i++)
	{
		printf(" ");
		slopes[i].print(i);
	}
	return 0;
}
