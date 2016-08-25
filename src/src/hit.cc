#include <cstring>
#include <cassert>
#include <cstdio>
#include <sstream>

#include "hit.h"
#include "config.h"

hit::hit(int32_t p)
{
	bam1_core_t::pos = p;
	strand = '.';
	xs = '.';
}

hit::hit(const hit &h)
	:bam1_core_t(h)
{
	rpos = h.rpos;
	qname = h.qname;
	strand = h.strand;
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
	/*
	if(bam_is_rev(b) == true) strand = '-';
	else strand = '+';
	*/

	strand = '.';
	if((flag & 0x10) <= 0 && (flag & 0x20) >= 1 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) strand = '+';		// F1R2
	if((flag & 0x10) >= 1 && (flag & 0x20) <= 0 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) strand = '-';		// R1F2
	if((flag & 0x10) <= 0 && (flag & 0x20) >= 1 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) strand = '-';		// F2R1
	if((flag & 0x10) >= 1 && (flag & 0x20) <= 0 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) strand = '+';		// R2F1

	if(strand_reverse == true && strand == '+') strand = '-';
	else if(strand_reverse == true && strand == '-') strand = '+';

	xs = '.';
	uint8_t *p = bam_aux_get(b, "XS");
	if(p && (*p) == 'A') xs = bam_aux2A(p);
	//if(xs == '-' || xs == '+') strand = xs;

	/*
	if((flag & 0x10) <= 0 && (flag & 0x20) >= 1 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) strand = '-';
	if((flag & 0x10) >= 1 && (flag & 0x20) <= 0 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) strand = '-';
	if((flag & 0x10) <= 0 && (flag & 0x20) >= 1 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) strand = '+';
	if((flag & 0x10) >= 1 && (flag & 0x20) <= 0 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) strand = '+';
	*/
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
	printf("Hit %s: [%d-%d), mpos = %d, cigar = %s, flag = %d, quality = %d, strand = %c, isize = %d\n", 
			qname.c_str(), pos, rpos, mpos, sstr.str().c_str(), flag, qual, strand, isize);

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

		if(k == 0 || k == n_cigar - 1) continue;
		if(bam_cigar_op(cigar[k]) != BAM_CREF_SKIP) continue;
		if(bam_cigar_op(cigar[k-1]) != BAM_CMATCH) continue;
		if(bam_cigar_op(cigar[k+1]) != BAM_CMATCH) continue;
		if(bam_cigar_oplen(cigar[k-1]) < min_flank_length) continue;
		if(bam_cigar_oplen(cigar[k+1]) < min_flank_length) continue;

		int32_t s = p - bam_cigar_oplen(cigar[k]);
		v.push_back(pack(s, p));
	}
    return 0;
}

int hit::get_mid_intervals(vector<int64_t> &vm, vector<int64_t> &vi, vector<int64_t> &vd) const
{
	vm.clear();
	vi.clear();
	vd.clear();
	int32_t p = pos;
    for(int k = 0; k < n_cigar; k++)
	{
		if (bam_cigar_type(bam_cigar_op(cigar[k]))&2)
		{
			p += bam_cigar_oplen(cigar[k]);
		}

		if(bam_cigar_op(cigar[k]) == BAM_CMATCH)
		{
			int32_t s = p - bam_cigar_oplen(cigar[k]);
			vm.push_back(pack(s, p));
		}

		if(bam_cigar_op(cigar[k]) == BAM_CINS)
		{
			vi.push_back(pack(p - 1, p + 1));
		}

		if(bam_cigar_op(cigar[k]) == BAM_CDEL)
		{
			int32_t s = p - bam_cigar_oplen(cigar[k]);
			vd.push_back(pack(s, p));
		}
	}
    return 0;
}

int hit::get_matched_intervals(vector<int64_t> &v) const
{
	vector<int64_t> vi, vd;
	return get_mid_intervals(v, vi, vd);
}

bool hit_compare_by_name(const hit &x, const hit &y)
{
	return x.qname < y.qname;
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
