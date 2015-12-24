#include <cstring>
#include <cassert>
#include <cstdio>
#include <sstream>

#include "hit.h"

hit::hit(bam1_t *b)
	:bam1_core_t(b->core)
{
	// fetch query name
	char buf[1024];
	memcpy(buf, bam_get_qname(b), l_qname);
	buf[l_qname] = '\0';
	qname = string(buf);

	// compute rpos
	rpos = pos + (int32_t)bam_cigar2rlen(n_cigar, bam_get_cigar(b));

	// copy cigar
	assert(n_cigar <= MAX_NUM_CIGAR);
	assert(n_cigar >= 1);
	memcpy(cigar, bam_get_cigar(b), 4 * n_cigar);

	// infer splice positions 
	n_spos = 0;
	infer_splice_positions();
}

hit::~hit()
{
}

int hit::print()
{
	// DEBUG
	if(n_spos <= 0) return 0;

	// get cigar string
	ostringstream sstr;
	for(int i = 0; i < n_cigar; i++)
	{
		sstr << bam_cigar_opchr(cigar[i]) << bam_cigar_oplen(cigar[i]);
		//printf("cigar %d: op = %c, length = %3d\n", i, bam_cigar_opchr(cigar[i]), bam_cigar_oplen(cigar[i]));
	}

	// print basic information
	printf("Hit %s: cigar = %s, flag = %d, quality = %d\n", qname.c_str(), sstr.str().c_str(), flag, qual);


	// print splice positions
	for(int i = 0; i < n_spos; i++)
	{
		printf(" splice %d: start = %7d, offset = %7d, splice = %7d\n", i, pos, spos[i] - pos, spos[i]);
	}
	return 0;
}

int hit::infer_splice_positions()
{
	int32_t p = pos;
    for(int k = 0; k < n_cigar; k++)
	{
		if (bam_cigar_type(bam_cigar_op(cigar[k]))&2)
			p += bam_cigar_oplen(cigar[k]);

		// must be flanked by matchings (with minimum length requirement of MIN_LEN_FLANK: TODO)
		if(k == 0 || k == n_cigar - 1) continue;
		if(bam_cigar_op(cigar[k]) != BAM_CREF_SKIP) continue;
		if(bam_cigar_op(cigar[k-1]) != BAM_CMATCH) continue;
		if(bam_cigar_op(cigar[k+1]) != BAM_CMATCH) continue;
		if(bam_cigar_oplen(cigar[k-1]) <= MIN_LEN_FLANK) continue;
		if(bam_cigar_oplen(cigar[k+1]) <= MIN_LEN_FLANK) continue;
		 
		assert(n_spos < MAX_NUM_CIGAR - 2);
		spos[n_spos++] = p - bam_cigar_oplen(cigar[k]);
		spos[n_spos++] = p;
	}
    return 0;
}

bool hit_compare_left(const hit &x, const hit &y)
{
	if(x.pos < y.pos) return true;
	else return false;
}

bool hit_compare_right(const hit &x, const hit &y)
{
	if(x.rpos < y.rpos) return true;
	return false;
}
