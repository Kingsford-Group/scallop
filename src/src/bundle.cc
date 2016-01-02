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
	infer_splice_boundaries();

	//sort(hits.begin(), hits.end(), hit_compare_left);
	check_left_ascending();
	count_prefix_left_hits();
	count_suffix_left_hits();

	/*
	sort(hits.begin(), hits.end(), hit_compare_right);
	check_right_ascending();

	count_prefix_right_hits();
	count_suffix_right_hits();

	infer_right_boundaries();

	sort(hits.begin(), hits.end(), hit_compare_left);
	check_left_ascending();
	*/

	infer_left_boundaries();

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

int bundle::count_suffix_left_hits()
{
	int p = 0;
	int q = 0;
	while(p < hits.size())
	{
		if(q >= hits.size() || hits[p].pos + hits_window_size < hits[q].pos)
		{
			hits[p].n_lhits += (int32_t)(q - p - 1);
			p++;
		}
		else
		{
			q++;
		}
	}
	return 0;
}

int bundle::count_prefix_left_hits()
{
	int p = hits.size() - 1;
	int q = hits.size() - 1;
	while(p >= 0)
	{
		if(q < 0 || hits[q].pos + hits_window_size < hits[p].pos)
		{
			hits[p].n_lhits += (int32_t)(p - q - 1);
			p--;
		}
		else
		{
			q--;
		}
	}
	return 0;
}


int bundle::count_suffix_right_hits()
{
	int p = 0;
	int q = 0;
	while(p < hits.size())
	{
		if(q >= hits.size() || hits[p].rpos + hits_window_size < hits[q].rpos)
		{
			hits[p].n_rhits += (int32_t)(q - p - 1);
			p++;
		}
		else
		{
			q++;
		}
	}
	return 0;
}

int bundle::count_prefix_right_hits()
{
	int p = hits.size() - 1;
	int q = hits.size() - 1;
	while(p >= 0)
	{
		if(q < 0 || hits[q].rpos + hits_window_size < hits[p].rpos)
		{
			hits[p].n_rhits += (int32_t)(p - q - 1);
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
				int type = p > 0 ? SPLICE_BOUNDARY_LEFT : SPLICE_BOUNDARY_RIGHT;
				int32_t pp = p > 0 ? p : 0 - p;
				boundary sp(type, pp, 1, hits[i].qual, hits[i].qual);
				m.insert(pair<int32_t, boundary>(p, sp));
			}
			else
			{
				m[p].count++;
				if(hits[i].qual < m[p].min_qual) m[p].min_qual = hits[i].qual;
				if(hits[i].qual > m[p].max_qual) m[p].max_qual = hits[i].qual;
			}
		}
	}

	map<int32_t, boundary>::iterator it;
	for(it = m.begin(); it != m.end(); it++)
	{
		if(it->second.count < min_splice_boundary_hits) continue;
		if(it->second.max_qual < min_max_splice_boundary_qual) continue;
		boundaries.push_back(it->second);
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

		binomial_distribution<> b(hits[i-1].n_lhits, 1.0/(2.0 * hits_window_size + 1.0));
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

		binomial_distribution<> b(hits[i-1].n_rhits, 1.0/(2.0 * hits_window_size + 1.0));
		double p = cdf(complement(b, sb.count - 1));

		sb.score = (uint32_t)(-10.0 * log10(p));

		if(sb.score < min_boundary_score) continue;

		boundaries.push_back(sb);
	}

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
