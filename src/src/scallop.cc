#include <cstdio>
#include <cassert>

#include "config.h"
#include "common.h"
#include "scallop.h"
#include "sam.h"

scallop::scallop()
{
}

scallop::~scallop()
{
}

int scallop::process(const char *bam_file)
{
	load(bam_file);
	solve();
	return 0;
}

int scallop::load(const char *bam_file)
{
    samFile *fn = sam_open(bam_file, "r");
    bam_hdr_t *h= sam_hdr_read(fn);
    bam1_t *b = bam_init1();

	bundle bd;
    while(sam_read1(fn, h, b) >= 0)
	{
		bam1_core_t &p = b->core;
		if((p.flag & 0x4) >= 1) continue;		// read is not mapped, TODO
		if((p.flag & 0x100) >= 1) continue;		// secondary alignment
		if(p.n_cigar < 1) continue;				// should never happen
		if(p.n_cigar > 7) continue;				// ignore hits with more than 7 cigar types
		//if(p.qual <= 4) continue;				// ignore hits with quality-score < 5
		if(bd.hits.size() > 0 && (bd.rpos + min_bundle_gap < p.pos || p.tid != bd.tid))
		{
			bundles.push_back(bd);
			bd.clear();

			// DEBUG
			// if(bundles.size() >= 10) break;
		}
		bd.add_hit(h, b);
    }

    bam_destroy1(b);
    bam_hdr_destroy(h);
    sam_close(fn);

	return 0;
}

int scallop::solve()
{
	for(int i = 0; i < bundles.size(); i++)
	{
		bundles[i].solve();

		// DEBUG
		if(bundles[i].chrm != "2L") continue;
		if(bundles[i].lpos < 870386) continue;
		if(bundles[i].rpos > 877183) continue;
		bundles[i].print();
	}
	return 0;
}
