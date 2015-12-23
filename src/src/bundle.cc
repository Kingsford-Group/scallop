#include <cassert>
#include <cstdio>

#include "kstring.h"
#include "bundle.h"

bundle::bundle()
{
	tid = -1;
	lpos = INT32_MAX;
	rpos = 0;
	n_hits = 0;
	m_hits = 1;
	hits = (bam1_t**)malloc(sizeof(bam1_t*) * m_hits);
	hits[0] = NULL;
}

bundle::~bundle()
{
	for(int i = 0; i < n_hits; i++) bam_destroy1(hits[i]);
	free(hits);
}

int bundle::add_hit(bam_hdr_t *h, bam1_t *b)
{
	// realloc memory if necessary
	assert(n_hits <= m_hits);
	if(n_hits == m_hits)
	{
		m_hits = n_hits * 2;
		hits = (bam1_t**)realloc(hits, m_hits * sizeof(bam1_t*)); 
	}

	// copy the hit
	hits[n_hits++] = bam_dup1(b);		// deep copy

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
		strcpy(chr, h->target_name[p.tid]);
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

int bundle::print()
{
	printf("tid = %4d, #hits = (%7d / %7d), range = %s:%d-%d\n", tid, n_hits, m_hits, chr, lpos, rpos);
	return 0;
}
