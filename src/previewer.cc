/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cassert>
#include <sstream>

#include "previewer.h"
#include "config.h"

previewer::previewer()
{
    sfn = sam_open(input_file.c_str(), "r");
    hdr = sam_hdr_read(sfn);
    b1t = bam_init1();
}

previewer::~previewer()
{
    bam_destroy1(b1t);
    bam_hdr_destroy(hdr);
    sam_close(sfn);
}

int previewer::preview()
{
	int total = 0;
	int first = 0;
	int second = 0;

    while(sam_read1(sfn, hdr, b1t) >= 0)
	{
		if(total >= max_preview_reads) break;

		bam1_core_t &p = b1t->core;

		if((p.flag & 0x4) >= 1) continue;										// read is not mapped
		if((p.flag & 0x100) >= 1 && use_second_alignment == false) continue;	// qstrandary alignment
		if(p.n_cigar > MAX_NUM_CIGAR) continue;									// ignore hits with more than 7 cigar types
		if(p.qual < min_mapping_quality) continue;								// ignore hits with small quality
		if(p.n_cigar < 1) continue;												// should never happen
		if((p.flag & 0x1) >= 1) continue;										// isoseq always single-end

		total++;

		hit ht(b1t);
		ht.set_tags(b1t);

		if(ht.xs == '.') continue;

		// predicted strand
		char xs = '.';
		if((ht.flag & 0x1) <= 0 && (ht.flag & 0x10) <= 0) xs = '+';
		if((ht.flag & 0x1) <= 0 && (ht.flag & 0x10) >= 1) xs = '-';

		if(xs == '.') continue;

		if(xs == ht.xs) first++;
		else second++;
	}

	vector<string> vv;
	vv.push_back("empty");
	vv.push_back("unstranded");
	vv.push_back("first");
	vv.push_back("second");

	int s1 = UNSTRANDED;
	int sp = first + second;
	if(sp >= min_preview_spliced_reads && first * 1.0 / sp > preview_infer_ratio) s1 = FR_FIRST;
	if(sp >= min_preview_spliced_reads && second * 1.0 / sp > preview_infer_ratio) s1 = FR_FIRST;

	if(verbose >= 1)
	{
		printf("preview: reads = %d, spliced reads = %d, first = %d, second = %d, inferred library_type = %s, given library_type = %s\n",
			total, sp, first, second, vv[s1 + 1].c_str(), vv[library_type + 1].c_str());
	}

	if(library_type == EMPTY) library_type = s1;

	return 0;
}
