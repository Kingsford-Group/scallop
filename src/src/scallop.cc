#include <cstdio>
#include <cassert>

#include "common.h"
#include "scallop.h"
#include "sam.h"

scallop::scallop(config *_conf)
	: conf(_conf)
{
}

scallop::~scallop()
{
}

int scallop::process(const char * bam_file)
{
    samFile *fn = sam_open(bam_file, "r");
    bam_hdr_t *h= sam_hdr_read(fn);
    bam1_t *b = bam_init1();

	bundle bd;
    while(sam_read1(fn, h, b) >= 0)
	{
		bam1_core_t &p = b->core;
		if((p.flag & 0x4) >= 1) continue;		// read is not mapped
		if((p.flag & 0x100) >= 1) continue;		// secondary alignment
		if(p.n_cigar < 1) continue;				// should never happen
		if(p.n_cigar > 7) continue;				// ignore hits with more than 7 cigar types
		if(bd.hits.size() > 0 && (bd.rpos + conf->min_bundle_gap < p.pos || p.tid != bd.tid))
		{
			bundles.push_back(bd);
			printf("bundle %8lu: ", bundles.size());
			bd.print();
			bd.clear();
		}
		bd.add_hit(h, b);
    }

    bam_destroy1(b);
    bam_hdr_destroy(h);
    sam_close(fn);

	return 0;
}
