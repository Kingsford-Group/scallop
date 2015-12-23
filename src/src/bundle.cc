#include <cassert>
#include <cstdio>

#include "kstring.h"
#include "bundle.h"

bundle::bundle()
{
	tid = -1;
	chr = "";
	lpos = INT32_MAX;
	rpos = 0;
}

bundle::~bundle()
{
}

int bundle::add_hit(bam_hdr_t *h, bam1_t *b)
{
	// store some information
	bam1_core_t &p = b->core;
	hits.push_back(p);

	// calcuate the boundaries on reference
	int32_t l = p.pos;
	int32_t r = p.pos + (int32_t)bam_cigar2rlen(p.n_cigar, bam_get_cigar(b));

	assert(r >= l);
	if(l < lpos) lpos = l;
	if(r > rpos) rpos = r;

	// set chromsome ID and name
	if(tid == -1)
	{
		char buf[1024];
		strcpy(buf, h->target_name[p.tid]);
		chr = string(buf);
		tid = p.tid;
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

int bundle::clear()
{
	hits.clear();
	tid = -1;
	chr = "";
	lpos = INT32_MAX;
	rpos = 0;
	return 0;
}

int bundle::print()
{
	printf("tid = %4d, #hits = %7ld, range = %s:%d-%d\n", tid, hits.size(), chr.c_str(), lpos, rpos);
	return 0;
}
