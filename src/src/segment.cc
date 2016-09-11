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
	select_slopes(0, bsets.size());
	build_partial_exons();
	return 0;
}

int segment::add_partial_exon(int k, const partial_exon &pe, double w)
{
	plist.push_back(k);
	pxs.push_back(pe);
	offsets.push_back(w);
	return 0;
}

int segment::build_bin_sets()
{
	bsets.clear();
	nbins.clear();
	for(int i = 0; i < pxs.size(); i++)
	{
		vector<int> bins;
		build_bins(pxs[i], bins, (int)(offsets[i]));
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
	for(SIMI it = lit; lit != imap->end(); it++)
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

	int n = 2;
	assert(bin_num % n == 0);
	int d = bin_num / n;

	for(int i = 0; i < bsets.size() - bin_num; i++)
	{
		vector<double> ave(n);
		vector<double> dev(n);
		vector<int> cnt(n);
		vector<int> score(n);
		vector<double> sigma(n);

		for(int k = 0; k < n; k++) 
		{
			compute_mean_dev(bsets, i + k * d, i + (k + 1) * d, ave[k], dev[k]);
			cnt[k] = ave[k] * d * slope_bin_size / average_read_length;
		}

		int32_t p1 = get_left_position(i + 1 * d);
		int32_t p2 = get_right_position(i + (n - 1) * d - 1);
		if(n <= 2) p2 = p1;

		for(int k = 0; k < n - 1; k++)
		{
			score[k] = compute_binomial_score(cnt[k] + cnt[k + 1], 0.5, cnt[k + 1]);
			sigma[k] = (ave[k + 1] - ave[k]) / dev[k];
		}
		score[n - 1] = compute_binomial_score(cnt[0] + cnt[n - 1], 0.5, cnt[n - 1]);
		sigma[n - 1] = (ave[n - 1] - ave[0]) / dev[0];

		if(score[n - 1] >= slope_min_score && sigma[n - 1] >= slope_min_sigma)
		{
			slope sp(SLOPE5END, i + d / 2, i + (n - 1) * d + d / 2, score[n - 1], sigma[n - 1]);
			sp.lpos = p1;
			sp.rpos = p2;
			seeds.push_back(sp);
		}

		for(int k = 0; k < n - 1; k++)
		{
			score[k] = compute_binomial_score(cnt[k] + cnt[k + 1], 0.5, cnt[k]);
			sigma[k] = (ave[k] - ave[k + 1]) / dev[k + 1];
		}
		score[n - 1] = compute_binomial_score(cnt[0] + cnt[n - 1], 0.5, cnt[0]);
		sigma[n - 1] = (ave[0] - ave[n - 1]) / dev[n - 1];

		if(score[n - 1] >= slope_min_score && sigma[n - 1] >= slope_min_sigma)
		{
			slope sp(SLOPE3END, i + d / 2, i + (n - 1) * d + d / 2, score[n - 1], sigma[n - 1]);
			sp.lpos = p1;
			sp.rpos = p2;
			seeds.push_back(sp);
		}
	}

	std::sort(seeds.begin(), seeds.end(), compare_slope_score);

	//for(int i = 0; i < seeds.size(); i++) seeds[i].print(i);

	return 0;
}

int segment::select_slopes(int si, int ti)
{
	int k = -1;
	for(int i = 0; i < seeds.size(); i++)
	{
		slope &x = seeds[i];
		if(x.lbin < si) continue;
		if(x.rbin > ti) continue;
		k = i;
		break;
	}

	if(k == -1) return 0;

	slope &x = seeds[k];

	//x.print(k);

	slopes.push_back(x);
	select_slopes(si, x.lbin);
	select_slopes(x.rbin, ti);

	return 0;
}

int segment::build_partial_exons()
{
	pexons.clear();
	if(pxs.size() == 0) return 0;
	sort(slopes.begin(), slopes.end(), compare_slope_pos);
	int i = 0, j = 0;
	while(i < pxs.size())
	{
		if(j == slopes.size())
		{
			pexons.push_back(pxs[i]);
			i++;
			continue;
		}

		slope &x = slopes[j];
		int32_t p = (x.type == SLOPE5END) ? x.lpos : x.rpos;
		int32_t p1 = pxs[i].lpos;
		int32_t p2 = pxs[i].rpos;

		//printf("i = %d, j = %d, p = %d, p1 = %d, p2 = %d, type = %d (%d, %d)\n", i, j, p, p1, p2, x.type, SLOPE5END, SLOPE3END);

		if(p == p1)
		{
			partial_exon pe = pxs[i];
			pe.ltype = START_BOUNDARY;
			pexons.push_back(pe);
			i++;
			j++;
		}
		else if(p == p2)
		{
			partial_exon pe = pxs[i];
			pe.rtype = END_BOUNDARY;
			pexons.push_back(pe);
			i++;
			j++;
		}
		else if(p < p1)
		{
			j++;
		}
		else if(p > p2)
		{
			pexons.push_back(pxs[i]);
			i++;
		}
		else if(pxs[i].ltype == START_BOUNDARY && x.type == SLOPE5END)
		{
			j++;
		}
		else if(pxs[i].rtype == END_BOUNDARY && x.type == SLOPE3END) 
		{
			j++;
		}
		else
		{
			partial_exon pe1(p1, p, pxs[i].ltype, x.type);
			partial_exon pe2(p, p2, x.type, pxs[i].rtype);
			evaluate_rectangle(*imap, pe1.lpos, pe1.rpos, pe1.ave, pe1.dev);
			evaluate_rectangle(*imap, pe2.lpos, pe2.rpos, pe2.ave, pe2.dev);
			pexons.push_back(pe1);
			pexons.push_back(pe2);
			i++;
			j++;
		}
	}
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
		return pxs[i].lpos + (x - s) * slope_bin_size;
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
		partial_exon &pe = pxs[i];
		if(x - s == nbins[i] - 1) return pe.rpos;
		else return pe.lpos + (x - s + 1) * slope_bin_size;
	}
	assert(false);
	return 0;
}

int segment::size() const
{
	return pxs.size();
}

int segment::print(int index) const
{
	printf("segment %d: \n", index);
	for(int i = 0; i < pxs.size(); i++)
	{
		printf("index = %d, offset = %.2lf, ", plist[i] + 1, offsets[i]);
		pxs[i].print(i);
	}
	for(int i = 0; i < seeds.size(); i++)
	{
		seeds[i].print(i);
	}
	for(int i = 0; i < slopes.size(); i++)
	{
		slopes[i].print(i);
	}
	/*
	for(int i = 0; i < pexons.size(); i++)
	{
		pexons[i].print(i);
	}
	*/
	return 0;
}
