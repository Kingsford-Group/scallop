#include <cassert>
#include <cstdio>
#include <map>

#include "kstring.h"
#include "bundle.h"

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

	infer_bridges();
	infer_left_boundaries();
	//infer_right_boundaries();
	add_start_boundary();
	add_end_boundary();

	remove_left_boundary_intervals();

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
	printf("tid = %d, #hits = %lu, range = %s:%d-%d\n", tid, hits.size(), chrm.c_str(), lpos, rpos);
	// print hits
	/*
	for(int i = 0; i < hits.size(); i++)
	{
		hits[i].print();
	}
	*/

	// print bridges 
	for(int i = 0; i < bridges.size(); i++)
	{
		bridges[i].print();
	}

	// print boundaries
	for(int i = 0; i < boundaries.size(); i++)
	{
		boundaries[i].print();
	}

	// print regions
	for(int i = 0; i < regions.size(); i++)
	{
		regions[i].print();
	}

	printf("\n");
	return 0;
}

int bundle::clear()
{
	tid = -1;
	chrm = "";
	lpos = INT32_MAX;
	rpos = 0;

	hits.clear();
	imap.clear();
	bridges.clear();
	boundaries.clear();
	regions.clear();

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

int bundle::remove_left_boundary_intervals()
{
	vector<int64_t> v;
	for(int i = 0; i < boundaries.size(); i++)
	{
		boundary &b = boundaries[i];
		if(b.type != LEFT_BOUNDARY) continue;

		int si = 0;
		int cnt = locate_hits(b.pos, si);
		if(cnt <= 0) continue;
		for(int j = si; j < cnt + si; j++)
		{
			if(hits[j].pos > hits[j].mpos) continue;
			if((hits[j].flag & 0x10) >= 1) continue;			

			hits[j].get_matched_intervals(v);
			for(int k = 0; k < v.size(); k++)
			{
				int32_t s = high32(v[k]);
				int32_t t = low32(v[k]);
				imap -= make_pair(ROI(s, t), 1);
			}
		}
	}
	return 0;
}

int bundle::locate_hits(int32_t p, int &li)
{
	li = -1;
	hit h(p);
	vector<hit>::iterator low = lower_bound(hits.begin(), hits.end(), h);
	vector<hit>::iterator up = upper_bound(hits.begin(), hits.end(), h);
	if(low == hits.end()) return 0;
	li = low - hits.begin();
	return (up - low);
}

int bundle::infer_bridges()
{
	map<int64_t, bridge> m;
	vector<int64_t> v;
	for(int i = 0; i < hits.size(); i++)
	{
		hits[i].get_splice_positions(v);
		if(v.size() == 0) continue;
		for(int k = 0; k < v.size(); k++)
		{
			int64_t p = v[k];
			if(m.find(p) == m.end()) 
			{
				bridge sp(SPLICE, p, 1, hits[i].qual, hits[i].qual);
				m.insert(pair<int64_t, bridge>(p, sp));
			}
			else
			{
				m[p].count++;
				if(hits[i].qual < m[p].min_qual) m[p].min_qual = hits[i].qual;
				if(hits[i].qual > m[p].max_qual) m[p].max_qual = hits[i].qual;
			}
		}
	}

	map<int64_t, bridge>::iterator it;
	for(it = m.begin(); it != m.end(); it++)
	{
		if(it->second.count < min_splice_boundary_hits) continue;
		if(it->second.max_qual < min_max_splice_boundary_qual) continue;
		bridges.push_back(it->second);
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

		int r = compute_overlap(imap, tpre);

		sb.score = compute_binomial_score(r, 1.0 / average_read_length, sb.count);

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

		int r = compute_overlap(imap, tpre);

		sb.score = compute_binomial_score(r, 1.0 / average_read_length, sb.count);

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
	set<int32_t> s;
	for(int i = 0; i < bridges.size(); i++)
	{
		s.insert(bridges[i].lpos);
		s.insert(bridges[i].rpos);
	}

	for(int i = 0; i < boundaries.size(); i++)
	{
		s.insert(boundaries[i].pos);
	}

	vector<int32_t> v(s.begin(), s.end());
	sort(v.begin(), v.end());

	if(v.size() <= 1) return 0;

	for(int i = 0; i < v.size() - 1; i++)
	{
		region r(v[i], v[i + 1], &imap);
		regions.push_back(r);
	}
	return 0;
}
