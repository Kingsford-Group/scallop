/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "region.h"
#include "config.h"
#include "util.h"
#include <algorithm>

using namespace std;

region::region(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype, const split_interval_map *_mmap, const split_interval_map *_imap)
	:lpos(_lpos), rpos(_rpos), mmap(_mmap), imap(_imap), ltype(_ltype), rtype(_rtype)
{

	build_join_interval_map();
	smooth_join_interval_map();
	//split_join_interval_map();
	build_partial_exons();
}

region::~region()
{}

int region::build_join_interval_map()
{
	jmap.clear();

	PSIMI pei = locate_boundary_iterators(*mmap, lpos, rpos);
	SIMI lit = pei.first, rit = pei.second;

	if(lit == mmap->end() || rit == mmap->end()) return 0;

	SIMI it = lit;
	while(true)
	{
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

int region::split_join_interval_map()
{
	if(lower(jmap.begin()->first) != lpos) return 0;
	if(upper(jmap.begin()->first) == rpos) return 0;

	PSIMI pei = locate_boundary_iterators(*mmap, lpos, rpos);
	SIMI lit = pei.first, rit = pei.second;
	if(lit == mmap->end() || rit == mmap->end()) return 0;

	int32_t min_split_middle_length = 10;
	int32_t min_split_boundary_length = 40;
	int32_t max_split_middle_coverage = 2;
	double min_split_boundary_coverage = 10;

	if(rpos - lpos < 100) return 0;

	int32_t p = lpos;
	for(SIMI it = lit; it != rit; it++)
	{
		int32_t p1 = lower(it->first);
		int32_t p2 = upper(it->first);
		int32_t cov = it->second;

		if(cov <= max_split_middle_coverage) continue;

		int32_t q = p1;
		bool b = true;
		if(q - p < min_split_middle_length) b = false;
		if(p - lpos < min_split_boundary_length) b = false;
		if(rpos - q < min_split_boundary_length) b = false;

		double ave1 = -1, ave2 = -1;
		double dev1 = -1, dev2 = -1;
		evaluate_rectangle(*mmap, lpos, p, ave1, dev1);
		evaluate_rectangle(*mmap, q, rpos, ave2, dev2);

		if(ave1 < min_split_boundary_coverage) b = false;
		if(ave2 < min_split_boundary_coverage) b = false;

		if(b == true)
		{
			printf("SPLIT: type = (%d, %d), pos = %d-%d, mid = %d-%d, ave = (%.3lf, %.3lf), dev = (%.3lf, %.3lf)\n", 
					ltype, rtype, lpos, rpos, p, q, ave1, ave2, dev1, dev2);

			jmap += make_pair(ROI(p, q), -1);
		}

		p = p2;
	}

	return 0;
}


int region::smooth_join_interval_map()
{
	int32_t gap = min_subregion_gap;
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

	PSIMI pei = locate_boundary_iterators(*mmap, p1, p2);
	SIMI it1 = pei.first, it2 = pei.second;
	if(it1 == mmap->end() || it2 == mmap->end()) return true;

	int32_t sum = compute_sum_overlap(*mmap, it1, it2);
	double ratio = sum * 1.0 / (p2 - p1);
	//printf(" region = [%d, %d), subregion [%d, %d), overlap = %.2lf\n", lpos, rpos, p1, p2, ratio);
	//if(ratio < min_subregion_overlap + max_intron_contamination_coverage) return true;
	if(ratio < min_subregion_overlap) return true;

	return false;
}

int region::build_partial_exons()
{
	pexons.clear();

	if(jmap.size() == 0) return 0;

	//printf("size = %lu, size2 = %lu, [%d, %d), [%d, %d)\n", jmap.size(), distance(jmap.begin(), jmap.end()), lower(jmap.begin()->first), upper(jmap.begin()->first), lpos, rpos);

	if(lower(jmap.begin()->first) == lpos && upper(jmap.begin()->first) == rpos)
	{
		partial_exon pe(lpos, rpos, ltype, rtype);
		evaluate_rectangle(*mmap, pe.lpos, pe.rpos, pe.ave, pe.dev);
		pexons.push_back(pe);
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
