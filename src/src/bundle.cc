#include <cassert>
#include <cstdio>
#include <map>

#include "kstring.h"
#include "bundle.h"

#include <boost/math/distributions/binomial.hpp>

using namespace boost::math;

bundle::bundle()
{
	tid = -1;
	chrm = "";
	lpos = INT32_MAX;
	rpos = 0;
}

bundle::~bundle()
{
}

int bundle::solve()
{
	// make sure all reads are sorted 
	check_left_ascending();

	build_interval_map();

	infer_splice_boundaries();
	infer_left_boundaries();
	//infer_right_boundaries();
	add_start_boundary();
	add_end_boundary();

	build_regions();

	return 0;
}

int bundle::add_hit(bam_hdr_t *h, bam1_t *b)
{
	// create and store new hit
	hit ht(b);
	hits.push_back(ht);

	// calcuate the boundaries on reference
	if(ht.pos < lpos) lpos = ht.pos;
	if(ht.rpos > rpos) rpos = ht.rpos;

	// set chromsome ID and name
	if(tid == -1)
	{
		tid = ht.tid;
		assert(chrm == "");
		char buf[1024];
		strcpy(buf, h->target_name[ht.tid]);
		chrm = string(buf);
	}
	assert(tid == ht.tid);

	return 0;
}

int bundle::print()
{
	printf("Bundle: ");
	printf("tid = %4d, #hits = %7lu, range = %s:%d-%d\n", tid, hits.size(), chrm.c_str(), lpos, rpos);
	// print hits
	/*
	for(int i = 0; i < hits.size(); i++)
	{
		hits[i].print();
	}
	*/

	// print boundaries
	for(int i = 0; i < boundaries.size(); i++)
	{
		boundaries[i].print();
	}

	printf("\n");
	return 0;
}

int bundle::clear()
{
	hits.clear();
	tid = -1;
	chrm = "";
	lpos = INT32_MAX;
	rpos = 0;
	return 0;
}

int bundle::build_interval_map()
{
	imap.clear();
	vector<int64_t> v;
	for(int i = 0; i < hits.size(); i++)
	{
		hits[i].get_matched_intervals(v);
		for(int k = 0; k < v.size(); k++)
		{
			int32_t s = high32(v[k]);
			int32_t t = low32(v[k]);
			imap += make_pair(ROI(s, t), 1);
		}
	}
	return 0;
}

int bundle::count_overlap_reads(int32_t p)
{
	imap_t::const_iterator it = imap.find(p);
	if(it == imap.end()) return 0;
	return it->second;
}

int bundle::infer_splice_boundaries()
{
	vector<int64_t> v;
	map<int64_t, uint32_t> m_count;
	map<int64_t, uint32_t> min_qual;
	map<int64_t, uint32_t> max_qual;

	for(int i = 0; i < hits.size(); i++)
	{
		hits[i].get_splice_positions(v);
		if(v.size() == 0) continue;
		for(int k = 0; k < v.size(); k++)
		{
			int64_t p = v[k];
			if(m_count.find(p) == m_count.end()) 
			{
				m_count.insert(pair<int64_t, uint32_t>(p, 1));
				max_qual.insert(pair<int64_t, uint32_t>(p, hits[i].qual));
				min_qual.insert(pair<int64_t, uint32_t>(p, hits[i].qual));
			}
			else
			{
				m_count[p]++;
				if(hits[i].qual > max_qual[p]) max_qual[p] = hits[i].qual;
				if(hits[i].qual < min_qual[p]) min_qual[p] = hits[i].qual;
			}
		}
	}

	map<int32_t, size_t> m;
	map<int64_t, uint32_t>::iterator it;
	for(it = m_count.begin(); it != m_count.end(); it++)
	{
		int64_t p = it->first;
		if(m_count[p] < min_splice_boundary_hits) continue;
		if(max_qual[p] < min_max_splice_boundary_qual) continue;

		int32_t s = high32(p);
		int32_t t = low32(p);
		size_t si = -1;
		size_t ti = -1;

		if(m.find(s) == m.end())
		{
			boundary b(SPLICE_BOUNDARY_LEFT, s, m_count[p], min_qual[p], max_qual[p]);
			si = boundaries.size();
			boundaries.push_back(b);
			m.insert(pair<int32_t, size_t>(s, si));
		}
		else
		{
			si = m[s];
			boundary &b = boundaries[si];
			assert(b.pos == s);
			b.count += m_count[p];
			if(min_qual[p] < b.min_qual) b.min_qual = min_qual[p];
			if(max_qual[p] > b.max_qual) b.max_qual = max_qual[p];
		}

		if(m.find(t) == m.end())
		{
			boundary b(SPLICE_BOUNDARY_RIGHT, t, m_count[p], min_qual[p], max_qual[p]);
			ti = boundaries.size();
			boundaries.push_back(b);
			m.insert(pair<int32_t, size_t>(t, ti));
		}
		else
		{
			ti = m[t];
			boundary &b = boundaries[ti];
			assert(b.pos == t);
			b.count += m_count[p];
			if(min_qual[p] < b.min_qual) b.min_qual = min_qual[p];
			if(max_qual[p] > b.max_qual) b.max_qual = max_qual[p];
		}

		assert(si != -1 && ti != -1);
		bridges.push_back(PT(si, ti));

	}
	return 0;
}

int bundle::infer_left_boundaries()
{
	int32_t pre = 0;
	int cnt = 1;
	for(int i = 1; i < hits.size(); i++)
	{
		if(hits[i].pos == pre)
		{
			cnt++;
			continue;
		}

		int32_t tpre = pre;
		int tcnt = cnt;
		
		cnt = 1;
		pre = hits[i].pos;

		// check the validity of the previous position
		if(tcnt < min_left_boundary_hits) continue;

		boundary sb(LEFT_BOUNDARY, tpre);
		for(int k = 0; k < tcnt; k++)
		{
			int j = i - 1 - k;
			if(hits[j].pos > hits[j].mpos) continue;
			if((hits[j].flag & 0x10) >= 1) continue;			
			sb.count++;
			if(hits[j].qual < sb.min_qual) sb.min_qual = hits[j].qual;
			if(hits[j].qual > sb.max_qual) sb.max_qual = hits[j].qual;
		}
		
		if(sb.count < min_left_boundary_hits) continue;
		if(sb.max_qual < min_max_left_boundary_qual) continue;

		int r = count_overlap_reads(tpre);

		binomial_distribution<> b(1.0 * r, 1.0 / average_read_length);
		double p = cdf(complement(b, sb.count - 1));

		sb.score = (uint32_t)(-10.0 * log10(p));

		if(sb.score < min_boundary_score) continue;

		boundaries.push_back(sb);
	}

	return 0;
}

int bundle::infer_right_boundaries()
{
	int32_t pre = 0;
	int cnt = 1;
	for(int i = 1; i < hits.size(); i++)
	{
		if(hits[i].rpos == pre)
		{
			cnt++;
			continue;
		}

		int32_t tpre = pre;
		int tcnt = cnt;
		
		cnt = 1;
		pre = hits[i].rpos;

		// check the validity of the previous position
		if(tcnt < min_right_boundary_hits) continue;

		boundary sb(RIGHT_BOUNDARY, tpre);
		for(int k = 0; k < tcnt; k++)
		{
			int j = i - 1 - k;
			if(hits[j].pos < hits[j].mpos) continue;
			if((hits[j].flag & 0x10) == 0) continue;			
			sb.count++;
			if(hits[j].qual < sb.min_qual) sb.min_qual = hits[j].qual;
			if(hits[j].qual > sb.max_qual) sb.max_qual = hits[j].qual;
		}
		
		if(sb.count < min_right_boundary_hits) continue;
		if(sb.max_qual < min_max_right_boundary_qual) continue;

		int32_t r = count_overlap_reads(tpre);

		binomial_distribution<> b(r, 1.0 / average_read_length);

		double p = cdf(complement(b, sb.count - 1));

		sb.score = (uint32_t)(-10.0 * log10(p));

		if(sb.score < min_boundary_score) continue;

		boundaries.push_back(sb);
	}

	return 0;
}

int bundle::add_start_boundary()
{
	for(int i = 0; i < boundaries.size(); i++)
	{
		boundary &b = boundaries[i];
		if(b.pos == lpos) return 0;
	}

	boundary b(START_BOUNDARY, lpos, 1, 100, 100);
	boundaries.push_back(b);

	return 0;
}

int bundle::add_end_boundary()
{
	for(int i = 0; i < boundaries.size(); i++)
	{
		boundary &b = boundaries[i];
		if(b.pos == rpos) return 0;
	}

	boundary b(END_BOUNDARY, rpos, 1, 100, 100);
	boundaries.push_back(b);

	return 0;
}

int bundle::check_left_ascending()
{
	for(int i = 1; i < hits.size(); i++)
	{
		int32_t p1 = hits[i - 1].pos;
		int32_t p2 = hits[i].pos;
		assert(p1 <= p2);
	}
	return 0;
}

int bundle::check_right_ascending()
{
	for(int i = 1; i < hits.size(); i++)
	{
		int32_t p1 = hits[i - 1].rpos;
		int32_t p2 = hits[i].rpos;
		assert(p1 <= p2);
	}
	return 0;
}

int bundle::build_regions()
{
	sort(boundaries.begin(), boundaries.end());
	return 0;
}
