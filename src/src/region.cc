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

	end5.clear();
	end3.clear();

	init();
	build_bins();
}

region::~region()
{}

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
	
	for(int i = 0; i < bins.size() - 1; i++)
	{
		bins[i] = bins[i] / bsize;
	}

	bins[bnum - 1] = bins[bnum - 1] / (rcore - lcore - bsize * bnum + bsize);

	return 0;
}

int region::add_boundary(int xi, int type)
{
	if(type == SLOPE5END) end5.push_back(xi);
	if(type == SLOPE3END) end3.push_back(xi);
	return 0;
}

int region::build_partial_exons(vector<partial_exon> &pexons)
{
	map<int, int> m;
	typedef pair<int, PI> PPI;
	for(int i = 0; i < end5.size(); i++)
	{
		//int p = lpos + end5[i] * slope_bin_size;
		int p = end5[i];
		if(m.find(p) == m.end()) m.insert(PI(p, START_BOUNDARY));
	}

	for(int i = 0; i < end3.size(); i++)
	{
		//int p = lpos + end3[i] * slope_bin_size;
		//if(end3[i] == bins.size()) p = rpos;
		int p = end3[i];
		if(m.find(p) == m.end()) m.insert(PI(p, END_BOUNDARY));
	}

	//if(m.find(lpos) == m.end()) m.insert(PI(lpos, ltype));
	//if(m.find(rpos) == m.end()) m.insert(PI(rpos, rtype));
	if(m.find(0) == m.end()) m.insert(PI(0, ltype));
	if(m.find(bins.size()) == m.end()) m.insert(PI(bins.size(), rtype));

	vector<PI> v(m.begin(), m.end());

	sort(v.begin(), v.end());

	vector<int> pp;
	for(int i = 0; i < v.size(); i++)
	{
		int p = lcore + v[i].first * slope_bin_size;
		if(v[i].first == bins.size()) p = rcore;
		assert(p >= lcore && p <= rcore);
		pp.push_back(p);
	}

	pexons.clear();
	for(int i = 0; i < v.size() - 1; i++)
	{
		partial_exon p(pp[i], pp[i + 1], v[i].second, v[i + 1].second);
		compute_mean_dev(bins, v[i].first, v[i + 1].first, p.ave, p.dev);
		pexons.push_back(p);
	}

	return 0;
}

int region::print(int index) const
{
	int32_t lc = compute_overlap(*imap, lcore);
	int32_t rc = compute_overlap(*imap, rcore - 1);
	printf("region %d: empty = %c, type = (%d, %d), pos = [%d, %d), core = [%d, %d), bins = %lu, coverage = (%d, %d)\n", 
			index, empty ? 'T' : 'F', ltype, rtype, lpos, rpos, lcore, rcore, bins.size(), lc, rc);
	return 0;
}
