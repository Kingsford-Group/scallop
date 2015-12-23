#include <cassert>
#include <cstdio>

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

int bundle::add_hit(bam_hdr_t *h, bam1_t *b)
{
	// create and store new hit
	hit ht(b);
	hits.push_back(ht);

	// calcuate the boundaries on reference
	bam1_core_t &p = b->core;
	int32_t l = p.pos;
	int32_t r = p.pos + (int32_t)bam_cigar2rlen(p.n_cigar, bam_get_cigar(b));

	assert(r >= l);
	if(l < lpos) lpos = l;
	if(r > rpos) rpos = r;

	// set chromsome ID and name
	if(tid == -1)
	{
		tid = p.tid;
		assert(chrm == "");
		char buf[1024];
		strcpy(buf, h->target_name[p.tid]);
		chrm = string(buf);
	}
	assert(tid == p.tid);

	// DEBUG
	/*
	if(chr == "X" && l >= 21512338 && r <= 21513626)
	{
		kstring_t ks = {0, 0, 0};
		sam_format1(h, b, &ks);
		printf("flag = %d, read: %s\n", (p.flag & 0x100), ks_str(&ks));
	}
	*/

	return 0;
}

int bundle::print()
{
	//printf("tid = %4d, #hits = %7lu, range = %s:%d-%d\n", tid, hits.size(), chrm.c_str(), lpos, rpos);
	for(int i = 0; i < hits.size(); i++)
	{
		//printf("hit %7d: ", i);
		hits[i].print();
	}

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

