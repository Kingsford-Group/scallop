/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __HIT_H__
#define __HIT_H__

#include <string>
#include <vector>

#include "htslib/sam.h"
#include "config.h"

using namespace std;


/*! @typedef
 @abstract Structure for core alignment information.
 @field  tid     chromosome ID, defined by bam_hdr_t
 @field  pos     0-based leftmost coordinate
 @field  bin     bin calculated by bam_reg2bin()
 @field  qual    mapping quality
 @field  l_qname length of the query name
 @field  flag    bitwise flag
 @field  n_cigar number of CIGAR operations
 @field  l_qseq  length of the query sequence (read)
 @field  mtid    chromosome ID of next read in template, defined by bam_hdr_t
 @field  mpos    0-based leftmost coordinate of next read in template

typedef struct {
    int32_t tid;
    int32_t pos;
    uint32_t bin:16, qual:8, l_qname:8;
    uint32_t flag:16, n_cigar:16;
    int32_t l_qseq;
    int32_t mtid;
    int32_t mpos;
    int32_t isize;
} bam1_core_t;
*/

class hit: public bam1_core_t
{
public:
	//hit(int32_t p);
	hit(bam1_t *b);
	hit(const hit &h);
	~hit();
	bool operator<(const hit &h) const;
	hit& operator=(const hit &h);

public:
	int32_t rpos;							// right position mapped to reference [pos, rpos)
	int32_t qlen;							// read length
	string qname;							// query name
	char strand;							// strandness
	char xs;								// XS aux in sam
	char ts;								// ts tag used in minimap2
	int32_t nh;								// NH aux in sam
	int32_t hi;								// HI aux in sam
	int32_t nm;								// NM aux in sam
	bool concordant;						// whether it is concordant
	vector<int64_t> spos;					// splice positions
	vector<int64_t> itvm;					// matched interval
	vector<int64_t> itvi;					// insert interval
	vector<int64_t> itvd;					// delete interval

public:
	int set_tags(bam1_t *b);
	int set_strand();
	int set_concordance();
	int print() const;
};

//inline bool hit_compare_by_name(const hit &x, const hit &y);

#endif
