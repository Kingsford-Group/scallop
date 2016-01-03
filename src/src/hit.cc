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
}

hit::~hit()
{
}

int hit::print()
{
	// get cigar string
	ostringstream sstr;
	for(int i = 0; i < n_cigar; i++)
	{
		sstr << bam_cigar_opchr(cigar[i]) << bam_cigar_oplen(cigar[i]);
		//printf("cigar %d: op = %c, length = %3d\n", i, bam_cigar_opchr(cigar[i]), bam_cigar_oplen(cigar[i]));
	}

	// print basic information
	printf("Hit %s: [%d-%d), cigar = %s, flag = %d, quality = %d\n", 
			qname.c_str(), pos, rpos, sstr.str().c_str(), flag, qual);

	return 0;
}

int hit::get_splice_positions(vector<int32_t> &v)
{
	v.clear();
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
		 
		v.push_back(p - bam_cigar_oplen(cigar[k]));
		v.push_back(0 - p);
	}
    return 0;
}

int hit::get_matched_intervals(vector<int64_t> & v)
{
	v.clear();

	int32_t p = pos;
    for(int k = 0; k < n_cigar; k++)
	{
		if (bam_cigar_type(bam_cigar_op(cigar[k]))&2)
			p += bam_cigar_oplen(cigar[k]);

		// must be flanked by matchings (with minimum length requirement of MIN_LEN_FLANK: TODO)
		if(bam_cigar_op(cigar[k]) != BAM_CMATCH) continue;

		int64_t t = p;
		int64_t s = p - bam_cigar_oplen(cigar[k]);

		int64_t x = (s << 32) | t;

		int32_t ss = (int32_t)(x >> 32);
		int32_t tt = (int32_t)((x << 32) >> 32);
		assert(ss == s);
		assert(tt == t);
		
		v.push_back(x);
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
