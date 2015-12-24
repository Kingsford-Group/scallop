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
	check();
	infer_splice_positions();
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

	// print splice positions
	for(int i = 0; i < splist.size(); i++)
	{
		splist[i]->print();
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

int bundle::infer_splice_positions()
{
	map<int32_t, splice_pos*> m;
	for(int i = 0; i < hits.size(); i++)
	{
		if(hits[i].n_spos <= 0) continue;
		for(int k = 0; k < hits[i].n_spos; k++)
		{
			int32_t p = hits[i].spos[k];
			if(m.find(p) == m.end()) 
			{
				splice_pos *sp = new splice_pos(p, 1, hits[i].qual, hits[i].qual);
				m.insert(pair<int32_t, splice_pos*>(p, sp));
			}
			else
			{
				assert(m[p]->pos == p);
				m[p]->count++;
				if(hits[i].qual < m[p]->min_qual) m[p]->min_qual = hits[i].qual;
				if(hits[i].qual > m[p]->max_qual) m[p]->max_qual = hits[i].qual;
			}
		}
	}

	map<int32_t, splice_pos*>::iterator it;
	for(it = m.begin(); it != m.end(); it++)
	{
		splist.push_back(it->second);
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
