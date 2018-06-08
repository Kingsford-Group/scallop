/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cstring>
#include <cassert>
#include <cstdio>
#include <sstream>
#include <cmath>

#include "hit.h"
#include "config.h"

/*
hit::hit(int32_t p)
{
	bam1_core_t::pos = p;
	strand = '.';
	xs = '.';
	ts = '.';
	hi = -1;
	nh = -1;
	nm = 0;
	qlen = 0;
	is_long_read = false;
	cigar = NULL;
}
*/

hit& hit::operator=(const hit &h)
{
	bam1_core_t::operator=(h);
	rpos = h.rpos;
	qlen = h.qlen;
	qname = h.qname;
	strand = h.strand;
	spos = h.spos;
	xs = h.xs;
	ts = h.ts;
	hi = h.hi;
	nm = h.nm;
	is_long_read = h.is_long_read;
	start_boundary = h.start_boundary;
	end_boundary = h.end_boundary;

	cigar = new uint32_t[h.n_cigar];
	memcpy(cigar, h.cigar, 4 * h.n_cigar);
	return *this;
}

hit::hit(const hit &h)
	:bam1_core_t(h)
{
	rpos = h.rpos;
	qlen = h.qlen;
	qname = h.qname;
	strand = h.strand;
	spos = h.spos;
	xs = h.xs;
	ts = h.ts;
	hi = h.hi;
	nm = h.nm;
	is_long_read = h.is_long_read;
	start_boundary = h.start_boundary;
	end_boundary = h.end_boundary;

	//printf("call copy constructor\n");
	cigar = new uint32_t[h.n_cigar];
	memcpy(cigar, h.cigar, 4 * h.n_cigar);
	//printf("n-cigar = %d, h.n_cigar = %d, size of cigar = %lu | %u\n", n_cigar, h.n_cigar, sizeof cigar, sizeof h.cigar);
}

hit::~hit()
{
	assert(cigar != NULL);
	delete[] cigar;
}

hit::hit(bam1_t *b)
	:bam1_core_t(b->core)
{
	// fetch query name
	char buf[1024];
	char *q = bam_get_qname(b);
	int l = strlen(q);
	memcpy(buf, q, l);
	buf[l] = '\0';
	qname = string(buf);

	/*
	string sub = qname.substr(0, 10);
	if(sub == "SRR1020625") is_long_read = false;
	else is_long_read = true;
	*/

	// compute rpos
	rpos = pos + (int32_t)bam_cigar2rlen(n_cigar, bam_get_cigar(b));
	qlen = (int32_t)bam_cigar2qlen(n_cigar, bam_get_cigar(b));

	// copy cigar
	assert(n_cigar <= MAX_NUM_CIGAR);
	assert(n_cigar >= 1);

	// allocate memery for cigar
	cigar = new uint32_t[n_cigar];
	memcpy(cigar, bam_get_cigar(b), 4 * n_cigar);

	//printf("call regular constructor\n");
}

int hit::set_tags(bam1_t *b)
{
	ts = '.';
	uint8_t *p0 = bam_aux_get(b, "ts");
	if(p0 && (*p0) == 'A') ts = bam_aux2A(p0);

	xs = '.';
	uint8_t *p1 = bam_aux_get(b, "XS");
	if(p1 && (*p1) == 'A') xs = bam_aux2A(p1);

	if(xs == '.' && ts != '.')
	{
		// convert ts to xs
		if((flag & 0x10) >= 1 && ts == '+') xs = '-';
		if((flag & 0x10) >= 1 && ts == '-') xs = '+';
		if((flag & 0x10) <= 0 && ts == '+') xs = '+';
		if((flag & 0x10) <= 0 && ts == '-') xs = '-';
	}

	hi = -1;
	uint8_t *p2 = bam_aux_get(b, "HI");
	if(p2 && (*p2) == 'C') hi = bam_aux2i(p2);

	nh = -1;
	uint8_t *p3 = bam_aux_get(b, "NH");
	if(p3 && (*p3) == 'C') nh = bam_aux2i(p3);

	nm = 0;
	uint8_t *p4 = bam_aux_get(b, "nM");
	if(p4 && (*p4) == 'C') nm = bam_aux2i(p4);

	uint8_t *p5 = bam_aux_get(b, "NM");
	if(p5 && (*p5) == 'C') nm = bam_aux2i(p5);

	return 0;
}

int hit::set_concordance()
{
	bool concordant = false;
	if((flag & 0x10) <= 0 && (flag & 0x20) >= 1 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) concordant = true;		// F1R2
	if((flag & 0x10) >= 1 && (flag & 0x20) <= 0 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) concordant = true;		// R1F2
	if((flag & 0x10) <= 0 && (flag & 0x20) >= 1 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) concordant = true;		// F2R1
	if((flag & 0x10) >= 1 && (flag & 0x20) <= 0 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) concordant = true;		// R2F1
	return 0;
}

int hit::set_strand()
{
	strand = xs;
	return 0;
}

int hit::set_ccsread_info(const ccsread_info &cci)
{
	start_boundary = end_boundary = false;
	if(strand == '+' && cci.fiveseen == true) start_boundary = true;
	if(strand == '+' && cci.threeseen == true) end_boundary = true;
	if(strand == '-' && cci.fiveseen == true) end_boundary = true;
	if(strand == '-' && cci.threeseen == true) start_boundary = true;
	return 0;
}

int hit::build_splice_positions()
{
	vector<PI32> vv;
	int32_t p = pos;
	int32_t q = 0;
	//uint8_t *seq = bam_get_seq(b);
    for(int k = 0; k < n_cigar; k++)
	{
		if (bam_cigar_type(bam_cigar_op(cigar[k]))&2)
			p += bam_cigar_oplen(cigar[k]);

		if (bam_cigar_type(bam_cigar_op(cigar[k]))&1)
			q += bam_cigar_oplen(cigar[k]);

		if(k == 0 || k == n_cigar - 1) continue;
		if(bam_cigar_op(cigar[k]) != BAM_CREF_SKIP) continue;

		/*
		if(bam_cigar_op(cigar[k-1]) != BAM_CMATCH) continue;
		if(bam_cigar_op(cigar[k+1]) != BAM_CMATCH) continue;
		if(bam_cigar_oplen(cigar[k-1]) < min_flank_length) continue;
		if(bam_cigar_oplen(cigar[k+1]) < min_flank_length) continue;
		*/

		int32_t s = p - bam_cigar_oplen(cigar[k]);
		vv.push_back(PI32(s, p));
	}

	spos.clear();
	if(vv.size() == 0) return 0;

	PI32 p1 = vv[0];
	for(int k = 1; k < vv.size(); k++)
	{
		if(vv[k].first == p1.second)
		{
			p1.second = vv[k].second;
		}
		else
		{
			spos.push_back(pack(p1.first, p1.second));
			p1 = vv[k];
		}
	}

	spos.push_back(pack(p1.first, p1.second));
	return 0;
}

bool hit::operator<(const hit &h) const
{
	if(qname < h.qname) return true;
	if(qname > h.qname) return false;
	if(hi != -1 && h.hi != -1 && hi < h.hi) return true;
	if(hi != -1 && h.hi != -1 && hi > h.hi) return false;
	return (pos < h.pos);
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
	printf("Hit %s: [%d-%d), mpos = %d, cigar = %s, flag = %d, quality = %d, strand = %c, xs = %c, ts = %c, isize = %d, qlen = %d, hi = %d, is-long = %c\n", 
			qname.c_str(), pos, rpos, mpos, sstr.str().c_str(), flag, qual, strand, xs, ts, isize, qlen, hi, is_long_read ? 'T' : 'F');

	printf(" start position (%d - )\n", pos);
	for(int i = 0; i < spos.size(); i++)
	{
		int64_t p = spos[i];
		int32_t p1 = high32(p);
		int32_t p2 = low32(p);
		printf(" splice position (%d - %d)\n", p1, p2);
	}
	printf(" end position (%d - )\n", rpos);

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
	v.clear();
	int32_t p1 = pos;
	int32_t p2 = -1;
	for(int i = 0; i < spos.size(); i++)
	{
		p2 = high32(spos[i]);
		v.push_back(pack(p1, p2));
		p1 = low32(spos[i]);
	}
	p2 = rpos;
	v.push_back(pack(p1, p2));
}

/*
inline bool hit_compare_by_name(const hit &x, const hit &y)
{
	if(x.qname < y.qname) return true;
	if(x.qname > y.qname) return false;
	return (x.pos < y.pos);
}

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

int hit::write_transcript(transcript &t)
{
	// assign attributes, and more
	t.seqname = "xxx";

	// add exons
	if(spos.size() == 0)
	{
		t.add_exon(pos, rpos);
	}
	else
	{
		int32_t l = pos;
		int64_t p = spos[0];
		int32_t p1 = low32(p);
		int32_t p2 = high32(p);
		l = p2;

		t.add_exon(l, p1);
		for(int k = 1; k < spos.size(); k++)
		{
			p = spos[k];
			p1 = low32(p);
			p2 = high32(p);
			t.add_exon(l, p1);
			l = p2;
		}

		t.add_exon(l, rpos);
	}

	// TODO
	t.strand = strand;
	return 0;
}
