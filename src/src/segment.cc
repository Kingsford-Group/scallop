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

		int32_t p2 = get_left_position(i + 1 * d);
		int32_t p3 = get_right_position(i + 2 * d - 1);

		int s1 = compute_binomial_score(xx + yy, 0.5, yy);
		int s2 = compute_binomial_score(yy + zz, 0.5, zz);
		int s3 = compute_binomial_score(xx + zz, 0.5, zz);

		double t1 = (ave2 - ave1) / dev1;
		double t2 = (ave3 - ave2) / dev2;
		double t3 = (ave3 - ave1) / dev1;

		bool b = (s3 >= slope_min_score);
		b = (b && (s1 >= 0.3 * slope_min_score));
		b = (b && (s2 >= 0.3 * slope_min_score));
		b = (b && (t3 >= slope_min_sigma));
		b = (b && (t1 >= slope_min_sigma * 0.5));
		b = (b && (t2 >= slope_min_sigma * 0.5));
	
		if(b == true)
		{
			slope sp(SLOPE5END, i + d, i + 2 * d, s3, t3);
			sp.lpos = p2;
			sp.rpos = p3;
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

		b = (s3 >= slope_min_score);
		b = (b && (s1 >= 0.3 * slope_min_score));
		b = (b && (s2 >= 0.3 * slope_min_score));
		b = (b && (t3 >= slope_min_sigma));
		b = (b && (t1 >= slope_min_sigma * 0.5));
		b = (b && (t2 >= slope_min_sigma * 0.5));

		if(b == true)
		{
			slope sp(SLOPE3END, i + d, i + 2 * d, s3, t3);
			sp.lpos = p2;
			sp.rpos = p3;
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

		if(p <= p1)
		{
			j++;
		}
		else if(p >= p2)
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
	for(int i = 0; i < pexons.size(); i++)
	{
		pexons[i].print(i);
	}
	return 0;
}
