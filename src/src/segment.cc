#include "segment.h"
#include "config.h"
#include "util.h"
#include "binomial.h"
#include <algorithm>

using namespace std;

segment::segment(const split_interval_map *_imap)
	:imap(_imap)
{
}

segment::~segment()
{}

int segment::clear()
{
	bsets.clear();
	nbins.clear();
	return 0;
}

int segment::build()
{
	build_bin_sets();
	build_seeds();
	return 0;
}

int segment::add_partial_exon(int k, const partial_exon &pe, double w)
{
	plist.push_back(k);
	pexons.push_back(pe);
	offsets.push_back(w);
	return 0;
}

int segment::build_bin_sets()
{
	bsets.clear();
	nbins.clear();
	for(int i = 0; i < pexons.size(); i++)
	{
		vector<int> bins;
		build_bins(pexons[i], bins, (int)(offsets[i]));
		bsets.insert(bsets.end(), bins.begin(), bins.end());
		nbins.push_back(bins.size());
	}
	return 0;
}

int segment::build_bins(const partial_exon &pe, vector<int> &bins, int offset)
{
	int32_t bsize = slope_bin_size;
	int bnum = ceil((pe.rpos - pe.lpos) * 1.0 / bsize);
	bins.resize(bnum, 0.0);

	SIMI lit, rit;
	tie(lit, rit) = locate_boundary_iterators(*imap, pe.lpos, pe.rpos);
	for(SIMI it = lit; ; it++)
	{
		int o = it->second;
		int32_t l = lower(it->first);
		int32_t r = upper(it->first);
		int bl = (l - pe.lpos) / bsize;
		int br = (r - 1 - pe.lpos) / bsize;

		assert(r > l);
		assert(br >= bl);
		assert(bl >= 0);
		assert(br < bins.size());

		for(int b = bl; b <= br; b++)
		{
			int x = pe.lpos + b * bsize;
			int y = pe.lpos + (b + 1) * bsize;
			if(b == bl) x = l;
			if(b == br) y = r;
			bins[b] += o * (y - x);
		}
		if(it == rit) break;
	}
	
	for(int i = 0; i < bins.size() - 1; i++)
	{
		bins[i] = bins[i] / bsize;
	}

	bins[bnum - 1] = bins[bnum - 1] / (pe.rpos - pe.lpos - bsize * bnum + bsize);

	for(int i = 0; i < bins.size(); i++) bins[i] += offset;

	return 0;
}

int segment::build_seeds()
{
	int bin_num = slope_bin_num;
	if(bsets.size() < bin_num) return 0;

	assert(bin_num % 3 == 0);
	int d = bin_num / 3;

	for(int i = 0; i < bsets.size() - bin_num; i++)
	{
		double ave1, ave2, ave3, dev1, dev2, dev3;
		compute_mean_dev(bsets, i + 0 * d, i + 1 * d, ave1, dev1);
		compute_mean_dev(bsets, i + 1 * d, i + 2 * d, ave2, dev2);
		compute_mean_dev(bsets, i + 2 * d, i + 3 * d, ave3, dev3);

		int xx = ave1 * d * slope_bin_size / average_read_length;
		int yy = ave2 * d * slope_bin_size / average_read_length;
		int zz = ave3 * d * slope_bin_size / average_read_length;

		int32_t p1 = get_left_position(i + 0 * d);
		int32_t p2 = get_left_position(i + 1 * d);
		int32_t p3 = get_left_position(i + 2 * d);
		int32_t p4 = get_left_position(i + 3 * d);

		int s1 = compute_binomial_score(xx + yy, 0.5, yy);
		int s2 = compute_binomial_score(yy + zz, 0.5, zz);
		int s3 = compute_binomial_score(xx + zz, 0.5, zz);

		double t1 = (ave2 - ave1) / dev1;
		double t2 = (ave3 - ave2) / dev2;
		double t3 = (ave3 - ave1) / dev1;

		if(s3 > slope_min_score && t3 >= slope_min_sigma)
		{
			slope sp(SLOPE5END, i + d, s3, t3, p2, p3);
			sp.score1 = s1;
			sp.score2 = s2;
			sp.sigma1 = t1;
			sp.sigma2 = t2;
			seeds.push_back(sp);
		}

		s1 = compute_binomial_score(xx + yy, 0.5, xx);
		s2 = compute_binomial_score(yy + zz, 0.5, yy);
		s3 = compute_binomial_score(xx + zz, 0.5, xx);

		t1 = (ave1 - ave2) / dev2;
		t2 = (ave2 - ave3) / dev3;
		t3 = (ave1 - ave3) / dev3;

		if(s3 > slope_min_score && t3 >= slope_min_sigma)
		{
			slope sp(SLOPE3END, i + d, s3, t3, p2, p3);
			sp.score1 = s1;
			sp.score2 = s2;
			sp.sigma1 = t1;
			sp.sigma2 = t2;
			seeds.push_back(sp);
		}
	}

	std::sort(seeds.begin(), seeds.end(), compare_slope_score);

	//for(int i = 0; i < seeds.size(); i++) seeds[i].print(i);

	return 0;
}

int segment::get_left_position(int x)
{
	int s = 0;
	for(int i = 0; i < nbins.size(); i++)
	{
		if(s + nbins[i] <= x) 
		{
			s += nbins[i];
			continue;
		}
		return pexons[i].lpos + (x - s) * slope_bin_size;
	}
	assert(false);
	return 0;
}

int segment::get_right_position(int x)
{
	int s = 0;
	for(int i = 0; i < nbins.size(); i++)
	{
		if(s + nbins[i] <= x)
		{
			s += nbins[i];
			continue;
		}
		partial_exon &pe = pexons[i];
		if(x - s == nbins[i] - 1) return pe.rpos;
		else return pe.lpos + (x - s + 1) * slope_bin_size;
	}
	assert(false);
	return 0;
}

int segment::print(int index) const
{
	printf("segment %d: \n", index);
	for(int i = 0; i < pexons.size(); i++)
	{
		printf("index = %d, offset = %.2lf, ", plist[i] + 1, offsets[i]);
		pexons[i].print(i);
	}
	for(int i = 0; i < seeds.size(); i++)
	{
		seeds[i].print(i);
	}
	return 0;
}
