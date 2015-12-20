#include <cstdio>
#include <cassert>

#include "common.h"
#include "scallop.h"
#include "sam.h"

scallop::scallop(config *_conf)
	: conf(_conf)
{}

scallop::~scallop()
{}

int scallop::process(const char * bam_file)
{
    samFile *fn = sam_open(bam_file, "r");
    bam_hdr_t *header = sam_hdr_read(fn);
    bam1_t *p = bam_init1();

	bundle b;
    while(sam_read1(fn, header, p) >= 0)
	{
		if(b.size() > 0 && b.rpos + conf->min_bundle_gap < p->core.pos)
		{
			bundles.push_back(b);
			printf("bundle %8ld: %8ld hits [%8d,%8d]\n", bundles.size(), b.size(), b.lpos, b.rpos);
			b.clear();
		}

		b.add_hit(p->core);
		continue;

		printf("template id = %d\n", p->core.tid);
		printf("position    = %d\n", p->core.pos);
		printf("bin         = %d\n", p->core.bin);
		printf("quality     = %d\n", p->core.qual);
		printf("qname len   = %d\n", p->core.l_qname);
		printf("qname       = %s\n", bam_get_qname(p));
		printf("flag        = %d\n", p->core.flag);
		printf("qseq len    = %d\n", p->core.l_qseq);
		printf("mate tid    = %d\n", p->core.mtid);
		printf("mate pos    = %d\n", p->core.mpos);
		printf("insert size = %d\n", p->core.isize);

		uint32_t *cigar = bam_get_cigar(p);
		printf("n_cigar     = %d\n", p->core.n_cigar);
		for(int k = 0; k < p->core.n_cigar; k++)
		{
			int32_t op = bam_cigar_op(cigar[k]);
			int32_t len = bam_cigar_oplen(cigar[k]);
			printf("cigar %1d     = (type = %d, length = %d)\n", k, op, len);
		}

		printf("\n");
    }

    bam_destroy1(p);
    bam_hdr_destroy(header);
    sam_close(fn);

	return 0;
}

//int scallop::process_line(bam1_t *p)
