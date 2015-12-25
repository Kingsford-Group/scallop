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
	check();
	count_prefix_hits();
	count_suffix_hits();
	infer_splice_boundaries();
	infer_start_boundaries();
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
	// print hits
	printf("tid = %4d, #hits = %7lu, range = %s:%d-%d\n", tid, hits.size(), chrm.c_str(), lpos, rpos);
	for(int i = 0; i < hits.size(); i++)
	{
		hits[i].print();
	}

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

int bundle::count_suffix_hits()
{
	int p = 0;
	int q = 0;
	while(p < hits.size())
	{
		if(q >= hits.size() || hits[p].pos + hits_window_size < hits[q].pos)
		{
			//printf("B: p = %d, window = %d, q = %d\n", hits[p].pos, hits_window_size, hits[q].pos);
			hits[p].n_hits += (int32_t)(q - p - 1);
			p++;
		}
		else
		{
			q++;
		}
	}
	return 0;
}

int bundle::count_prefix_hits()
{
	int p = hits.size() - 1;
	int q = hits.size() - 1;
	while(p >= 0)
	{
		if(q < 0 || hits[q].pos + hits_window_size < hits[p].pos)
		{
			hits[p].n_hits += (int32_t)(p - q - 1);
			p--;
		}
		else
		{
			q--;
		}
	}
	return 0;
}

int bundle::infer_splice_boundaries()
{
	map<int32_t, boundary> m;
	for(int i = 0; i < hits.size(); i++)
	{
		if(hits[i].n_spos <= 0) continue;
		for(int k = 0; k < hits[i].n_spos; k++)
		{
			int32_t p = hits[i].spos[k];
			if(m.find(p) == m.end()) 
			{
				boundary sp(SPLICE_BOUNDARY, p, 1, hits[i].qual, hits[i].qual);
				m.insert(pair<int32_t, boundary>(p, sp));
			}
			else
			{
				assert(m[p].pos == p);
				m[p].count++;
				if(hits[i].qual < m[p].min_qual) m[p].min_qual = hits[i].qual;
				if(hits[i].qual > m[p].max_qual) m[p].max_qual = hits[i].qual;
			}
		}
	}

	map<int32_t, boundary>::iterator it;
	for(it = m.begin(); it != m.end(); it++)
	{
		if(it->second.count < min_start_boundary_hits) continue;
		if(it->second.max_qual < min_max_start_boundary_qual) continue;
		boundaries.push_back(it->second);
	}
	return 0;
}

int bundle::infer_start_boundaries()
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
		if(tcnt < min_start_boundary_hits) continue;

		boundary sb(START_BOUNDARY, tpre);
		for(int k = 0; k < tcnt; k++)
		{
			int j = i - 1 - k;
			if(hits[j].pos > hits[j].mpos) continue;		// consider boundary of a SEGMENT, not a read
			sb.count++;
			if(hits[j].qual < sb.min_qual) sb.min_qual = hits[j].qual;
			if(hits[j].qual > sb.max_qual) sb.max_qual = hits[j].qual;
		}
		
		if(sb.count < min_start_boundary_hits) continue;
		if(sb.max_qual < min_max_start_boundary_qual) continue;

		binomial_distribution<> b(hits[i].n_hits, 1.0/(2.0 * hits_window_size + 1.0));
		double p = cdf(complement(b, sb.count - 1));

		sb.score = (uint32_t)(-10.0 * log10(p));

		if(sb.score < min_boundary_score) continue;

		boundaries.push_back(sb);
	}

	return 0;
}

int bundle::check()
{
	// guarantee that all reads are sorted
	for(int i = 1; i < hits.size(); i++)
	{
		int32_t p1 = hits[i - 1].pos;
		int32_t p2 = hits[i].pos;
		assert(p1 <= p2);
	}
	return 0;
}
