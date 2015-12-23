#include <stdio.h>
#include <assert.h>
#include <malloc.h>

#include "common.h"
#include "scallop.h"
#include "sam.h"

scallop::scallop(config *_conf)
	: conf(_conf)
{
	n_bundles = 0;
	m_bundles = 1;
	bundles = (bundle**)malloc(sizeof(bundle*) * m_bundles);
}

scallop::~scallop()
{
	for(int i = 0; i < n_bundles; i++) delete bundles[i];
	free(bundles);
}

int scallop::add_bundle(bundle *bd)
{
	assert(n_bundles <= m_bundles);
	if(n_bundles == m_bundles) 
	{
		m_bundles = n_bundles * 2;
		bundles = (bundle**)realloc(bundles, m_bundles * sizeof(bundle*));
	}

	bundles[n_bundles++] = bd;	// shallow copy!
	return 0;
}

int scallop::process(const char * bam_file)
{
    samFile *fn = sam_open(bam_file, "r");
    bam_hdr_t *h= sam_hdr_read(fn);
    bam1_t *b = bam_init1();

	bundle *bd = new bundle;
    while(sam_read1(fn, h, b) >= 0)
	{
		bam1_core_t &p = b->core;
		if((p.flag & 0x4) >= 1) continue;		// read is not mapped
		if((p.flag & 0x100) >= 1) continue;		// secondary alignment
		if(bd->n_hits > 0 && (bd->rpos + conf->min_bundle_gap < p.pos || p.tid != bd->tid))
		{
			add_bundle(bd);
			printf("bundle %8d: ", n_bundles);
			bd->print();
			bd = new bundle;
		}
		bd->add_hit(h, b);
    }

	delete bd;
    bam_destroy1(b);
    bam_hdr_destroy(h);
    sam_close(fn);

	return 0;
}

//int scallop::process_line(bam1_t *p)
