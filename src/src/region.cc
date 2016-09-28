#include "region.h"
#include "config.h"
#include "util.h"
#include "binomial.h"
#include <algorithm>

using namespace std;

region::region(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype, const split_interval_map *_mmap, const split_interval_map *_imap)
	:lpos(_lpos), rpos(_rpos), mmap(_mmap), imap(_imap), ltype(_ltype), rtype(_rtype)
{

	build_join_interval_map();
	smooth_join_interval_map();
	build_partial_exons();
}

region::~region()
{}

int region::build_join_interval_map()
{
	jmap.clear();

	SIMI lit, rit;
	tie(lit, rit) = locate_boundary_iterators(*mmap, lpos, rpos);
	if(lit == mmap->end() || rit == mmap->end()) return 0;

	SIMI it = lit;
	while(true)
	{
		// TODO
		//if(it->second >= 2) 
		jmap += make_pair(it->first, 1);
		if(it == rit) break;
		it++;
	}

	for(JIMI it = jmap.begin(); it != jmap.end(); it++)
	{
		assert(it->second == 1);
	}

	return 0;
}

int region::smooth_join_interval_map()
{
	int32_t gap = min_subregion_gap;

	/*
	bool b1 = false, b2 = false;
	if(ltype == START_BOUNDARY) b1 = true;
	if(ltype == RIGHT_SPLICE) b1 = true;
	if(rtype == END_BOUNDARY) b2 = true;
	if(rtype == LEFT_SPLICE) b2 = true;
	if(b1 == true && b2 == true) gap = 2 * min_subregion_gap;
	*/

	vector<PI32> v;
	int32_t p = lpos;
	for(JIMI it = jmap.begin(); it != jmap.end(); it++)
	{
		int32_t p1 = lower(it->first);
		int32_t p2 = upper(it->first);
		assert(p1 >= p);
		assert(p2 > p1);
		if(p1 - p <= gap) v.push_back(PI32(p, p1));
		p = p2;
	}

	if(p < rpos && rpos - p <= gap) v.push_back(PI32(p, rpos));

	for(int i = 0; i < v.size(); i++)
	{
		jmap += make_pair(ROI(v[i].first, v[i].second), 1);
	}

	for(JIMI it = jmap.begin(); it != jmap.end(); it++)
	{
		assert(it->second == 1);
	}

	return 0;
}

bool region::empty_subregion(int32_t p1, int32_t p2)
{
	assert(p1 < p2);
	assert(p1 >= lpos && p2 <= rpos);

	//printf(" region = [%d, %d), subregion [%d, %d), length = %d\n", lpos, rpos, p1, p2, p2 - p1);
	if(p2 - p1 < min_subregion_length) return true;

	SIMI it1, it2;
	tie(it1, it2) = locate_boundary_iterators(*mmap, p1, p2);
	if(it1 == mmap->end() || it2 == mmap->end()) return true;

	int ladders = distance(it1, it2);
	if(ladders < min_subregion_ladders) return false;

	int32_t sum = compute_sum_overlap(*mmap, it1, it2);
	double ratio = sum * 1.0 / (p2 - p1);
	//printf(" region = [%d, %d), subregion [%d, %d), overlap = %.2lf\n", lpos, rpos, p1, p2, ratio);
	if(ratio < min_subregion_overlap) return true;

	int32_t indel = 0;
	SIMI jit1 = imap->lower_bound(ROI(p1, p1 + 1));
	SIMI jit2 = imap->upper_bound(ROI(p2 - 1, p2));
	for(SIMI jit = jit1; jit != jit2; jit++) indel += jit->second;

	//printf(" region = [%d, %d), subregion [%d, %d), indel = %d\n", lpos, rpos, p1, p2, indel);
	if(indel * 1.0 / (p2 - p1) > max_indel_ratio) return true;

	return false;
}

int32_t region::identify_boundary(bool tag)
{
	if(identify_extra_boundary == false) return -1;
	if(lower(jmap.begin()->first) != lpos || upper(jmap.begin()->first) != rpos) return -1;

	typedef pair<int32_t, int64_t> PI64;
	vector<PI64> v;
	SIMI lit, rit;
	tie(lit, rit) = locate_boundary_iterators(*mmap, lpos, rpos);
	if(lit == mmap->end() || rit == mmap->end()) return -1;

	v.push_back(PPI(lpos, 0));
	for(SIMI it = lit; it != rit; it++)
	{
		int32_t s = lower(it->first);
		int32_t t = upper(it->first);
		int k = v.size() - 1;
		int64_t x = v[k].second + (t - s) * it->second;
		v.push_back(PI64(t, x));
	}

	int64_t sum = v[v.size() - 1].second;

	int32_t pos = -1;
	uint32_t score = 0;
	double ave1 = 0, ave2 = 0;
	for(int i = 0; i < v.size(); i++)
	{
		int32_t p = v[i].first;
		int32_t len1 = p - lpos;
		int32_t len2 = rpos - p;
		
		if(len1 < min_boundary_length) continue;
		if(len2 < min_boundary_length) continue;

		int n1 = v[i].second / average_read_length + 5;
		int n2 = (sum - v[i].second) / average_read_length + 5;

		uint32_t s = 0;
		if(tag == true)
		{
			double pr = len1 * 1.0 / (len1 + len2);
			s = compute_binomial_score(n1 + n2, pr, n1);
		}
		else
		{
			double pr = len2 * 1.0 / (len1 + len2);
			s = compute_binomial_score(n1 + n2, pr, n2);
		}

		if(score > s) continue;

		score = s;
		pos = p;
		ave1 = v[i].second * 1.0 / len1;
		ave2 = (sum - v[i].second) * 1.0 / len2;
	}

	if(pos == -1) return -1;

	int32_t p = lpos;
	double var1 = 0, var2 = 0;
	for(SIMI it = lit; it != rit; it++)
	{
		int32_t s = lower(it->first);
		int32_t t = upper(it->first);
		assert(s >= p);

		if(t <= pos)
		{
			var1 += (s - p) * ave1 * ave1;
			var1 += (t - s) * (it->second - ave1) * (it->second - ave1);
		}
		else
		{
			var2 += (s - p) * ave2 * ave2;
			var2 += (t - s) * (it->second - ave2) * (it->second - ave2);
		}
		p = t;
	}
	var2 += (rpos - p) * ave2 * ave2;

	int len1 = pos - lpos;
	int len2 = rpos - pos;
	double dev1 = sqrt(var1 / len1);
	double dev2 = sqrt(var2 / len2);

	if(score < min_boundary_score) return -1;
	if(tag == true && ave1 < ave2 + min_boundary_sigma * dev2) return -1;
	if(tag == false && ave2 < ave1 + min_boundary_sigma * dev1) return -1;

	printf("new boundary: %s-position: region = %d-%d, score = %d, pos = %d, len = (%d, %d) ave = (%.1lf, %.1lf), dev = (%.1lf, %.1lf)\n", 
			tag ? "end" : "start", lpos, rpos, score, pos, len1, len2, ave1, ave2, dev1, dev2);

	return pos;
}

int region::build_partial_exons()
{
	pexons.clear();

	if(jmap.size() == 0) return 0;

	//printf("size = %lu, size2 = %lu, [%d, %d), [%d, %d)\n", jmap.size(), distance(jmap.begin(), jmap.end()), lower(jmap.begin()->first), upper(jmap.begin()->first), lpos, rpos);

	if(lower(jmap.begin()->first) == lpos && upper(jmap.begin()->first) == rpos)
	{
		if(ltype == START_BOUNDARY || rtype == END_BOUNDARY)
		{
			partial_exon pe(lpos, rpos, ltype, rtype);
			evaluate_rectangle(*mmap, pe.lpos, pe.rpos, pe.ave, pe.dev);
			pexons.push_back(pe);

			//printf("whole region 1: "); pe.print(9);

			return 0;
		}

		int32_t p1 = identify_boundary(true);
		int32_t p2 = identify_boundary(false);

		if(p1 < 0 && p2 < 0)
		{
			partial_exon pe(lpos, rpos, ltype, rtype);
			evaluate_rectangle(*mmap, pe.lpos, pe.rpos, pe.ave, pe.dev);
			pexons.push_back(pe);

			//printf("whole region 2: "); pe.print(9);

			return 0;
		}
		else if(p1 > 0)
		{
			partial_exon pe1(lpos, p1, ltype, END_BOUNDARY);
			partial_exon pe2(p1, rpos, MIDDLE_CUT, rtype);
			evaluate_rectangle(*mmap, pe1.lpos, pe1.rpos, pe1.ave, pe1.dev);
			evaluate_rectangle(*mmap, pe2.lpos, pe2.rpos, pe2.ave, pe2.dev);
			pexons.push_back(pe1);
			pexons.push_back(pe2);
		}
		else if(p2 > 0)
		{
			partial_exon pe1(lpos, p2, ltype, MIDDLE_CUT);
			partial_exon pe2(p2, rpos, START_BOUNDARY, rtype);
			evaluate_rectangle(*mmap, pe1.lpos, pe1.rpos, pe1.ave, pe1.dev);
			evaluate_rectangle(*mmap, pe2.lpos, pe2.rpos, pe2.ave, pe2.dev);
			pexons.push_back(pe1);
			pexons.push_back(pe2);
		}
		return 0;
	}

	if(ltype == RIGHT_SPLICE && jmap.find(ROI(lpos, lpos + 1)) == jmap.end())
	{
		partial_exon pe(lpos, lpos + 1, ltype, END_BOUNDARY);
		pe.ave = 1.0;
		pe.dev = 1.0;
		pexons.push_back(pe);
	}

	for(JIMI it = jmap.begin(); it != jmap.end(); it++)
	{
		int32_t p1 = lower(it->first);
		int32_t p2 = upper(it->first);
		assert(p1 < p2);
		
		bool b = empty_subregion(p1, p2);

		//printf(" subregion [%d, %d), empty = %c\n", p1, p2, b ? 'T' : 'F');

		if(p1 == lpos && ltype == RIGHT_SPLICE) b = false;
		if(p2 == rpos && rtype == LEFT_SPLICE) b = false;

		if(b == true) continue;

		int lt = (p1 == lpos) ? ltype : START_BOUNDARY;
		int rt = (p2 == rpos) ? rtype : END_BOUNDARY;

		partial_exon pe(p1, p2, lt, rt);
		evaluate_rectangle(*mmap, pe.lpos, pe.rpos, pe.ave, pe.dev);
		pexons.push_back(pe);
	}

	if(rtype == LEFT_SPLICE && jmap.find(ROI(rpos - 1, rpos)) == jmap.end())
	{
		partial_exon pe(rpos - 1, rpos, START_BOUNDARY, rtype);
		pe.ave = 1.0;
		pe.dev = 1.0;
		pexons.push_back(pe);
	}

	return 0;
}

bool region::left_inclusive()
{
	if(pexons.size() == 0) return false;
	if(pexons[0].lpos == lpos) return true;
	else return false;
}

bool region::right_inclusive()
{
	if(pexons.size() == 0) return false;
	if(pexons[pexons.size() - 1].rpos == rpos) return true;
	else return false;
}

int region::print(int index) const
{
	int32_t lc = compute_overlap(*mmap, lpos);
	int32_t rc = compute_overlap(*mmap, rpos - 1);
	printf("region %d: partial-exons = %lu, type = (%d, %d), pos = [%d, %d), boundary coverage = (%d, %d)\n", 
			index, pexons.size(), ltype, rtype, lpos, rpos, lc, rc);

	/*
	for(JIMI it = jmap.begin(); it != jmap.end(); it++)
	{
		printf(" [%d, %d) -> %d\n", lower(it->first), upper(it->first), it->second);
	}
	*/

	return 0;
}
