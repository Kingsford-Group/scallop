#include <cstring>
#include <cassert>
#include <cstdio>
#include <sstream>

#include "hit.h"

hit::hit(int32_t p)
{
	bam1_core_t::pos = p;
	xs = '.';
}

hit::hit(const hit &h)
	:bam1_core_t(h)
{
	rpos = h.rpos;
	qname = h.qname;
	xs = h.xs;
	for(int i = 0; i < MAX_NUM_CIGAR; i++)
	{
		cigar[i] = h.cigar[i];
	}
}

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

	// get strandness
	uint8_t *p = bam_aux_get(b, "XS");
	if(p && (*p) == 'A')
	{
		xs = bam_aux2A(p);
	}
}

hit::~hit()
{
}

bool hit::operator<(const hit &h) const
{
	if(pos < h.pos) return true;
	else return false;
}

int hit::print() const
{
	// get cigar string
	ostringstream sstr;
	for(int i = 0; i < n_cigar; i++)
	{
		sstr << bam_cigar_opchr(cigar[i]) << bam_cigar_oplen(cigar[i]);
		//printf("cigar %d: op = %c, length = %3d\n", i, bam_cigar_opchr(cigar[i]), bam_cigar_oplen(cigar[i]));
	}

	// print basic information
	printf("Hit %s: [%d-%d), mpos = %d, cigar = %s, flag = %d, quality = %d, strand = %c\n", 
			qname.c_str(), pos, rpos, mpos, sstr.str().c_str(), flag, qual, xs);

	return 0;
}

int hit::get_splice_positions(vector<int64_t> &v) const
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

		int32_t s = p - bam_cigar_oplen(cigar[k]);
		v.push_back(pack(s, p));
	}
    return 0;
}

int hit::get_matched_intervals(vector<int64_t> & v) const
{
	v.clear();
	int32_t p = pos;
    for(int k = 0; k < n_cigar; k++)
	{
		if (bam_cigar_type(bam_cigar_op(cigar[k]))&2)
			p += bam_cigar_oplen(cigar[k]);

		// must be flanked by matchings (with minimum length requirement of MIN_LEN_FLANK: TODO)
		if(bam_cigar_op(cigar[k]) != BAM_CMATCH) continue;

		int32_t s = p - bam_cigar_oplen(cigar[k]);
		v.push_back(pack(s, p));
	}

    return 0;
}

/*
inline bool hit_compare_left(const hit &x, const hit &y)
{
	if(x.pos < y.pos) return true;
	else return false;
}

inline bool hit_compare_right(const hit &x, const hit &y)
{
	if(x.rpos < y.rpos) return true;
	return false;
}
*/
