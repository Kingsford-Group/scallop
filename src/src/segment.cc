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

	assert(bin_num % 2 == 0);
	int d = bin_num / 2;

	for(int i = 0; i < bsets.size() - bin_num; i++)
	{
		double ave1, ave2, dev1, dev2;
		compute_mean_dev(bsets, i, i + d, ave1, dev1);
		compute_mean_dev(bsets, i + d, i + 2 * d, ave2, dev2);
		int xx = ave1 * d * slope_bin_size / average_read_length;
		int yy = ave2 * d * slope_bin_size / average_read_length;
		int score5 = compute_binomial_score(xx + yy, 0.5, yy);
		int score3 = compute_binomial_score(xx + yy, 0.5, xx);

		if(score5 > slope_min_score)
		{
			slope s5(SLOPE5END, i + d, score5, (ave2 - ave1) / dev1);
			seeds.push_back(s5);
		}

		if(score3 > slope_min_score)
		{
			slope s3(SLOPE3END, i + d, score3, (ave1 - ave2) / dev2);
			seeds.push_back(s3);
		}
	}

	std::sort(seeds.begin(), seeds.end(), compare_slope_score);

	//for(int i = 0; i < seeds.size(); i++) seeds[i].print(i);

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
