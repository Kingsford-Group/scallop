#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "sam.h"

int test(const char * sam)
{
    samFile *in = sam_open(sam, "r");
    bam_hdr_t *header = sam_hdr_read(in);
    bam1_t *p = bam_init1();

    while(sam_read1(in, header, p) >= 0)
	{
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
			uint32_t op = bam_cigar_op(cigar[k]);
			uint32_t len = bam_cigar_oplen(cigar[k]);
			printf("cigar %1d     = (type = %d, length = %d)\n", k, op, len);
		}

		printf("\n");
    }

    bam_destroy1(p);
    bam_hdr_destroy(header);
    sam_close(in);

    return 0;
}

