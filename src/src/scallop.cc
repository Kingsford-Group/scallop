#include <cstdio>
#include <cassert>

#include "config.h"
#include "scallop.h"
#include "assembler.h"

scallop::scallop()
{
}

scallop::~scallop()
{
}

int scallop::process(const string &file)
{
	string s = file.substr(file.size() - 3, 3);
	if(s == "bam" || s == "sam") solve_bam(file);
	else if(s == "gtf") solve_gtf(file);
	else solve_example(file);
	return 0;
}

int scallop::solve_bam(const string &file)
{
    samFile *fn = sam_open(file.c_str(), "r");
    bam_hdr_t *h= sam_hdr_read(fn);
    bam1_t *b = bam_init1();

	ofstream fout1("greedy.gtf");
	ofstream fout2("iterat.gtf");

	int count = 0;
	bundle_base bb;
    while(sam_read1(fn, h, b) >= 0)
	{
		bam1_core_t &p = b->core;
		if((p.flag & 0x4) >= 1) continue;		// read is not mapped, TODO
		if((p.flag & 0x100) >= 1) continue;		// secondary alignment
		if(p.n_cigar < 1) continue;				// should never happen
		if(p.n_cigar > 7) continue;				// ignore hits with more than 7 cigar types
		//if(p.qual <= 4) continue;				// ignore hits with quality-score < 5
		if(bb.get_num_hits() > 0 && (bb.get_rpos() + min_bundle_gap < p.pos || p.tid != bb.get_tid()))
		{
			if(bb.get_num_hits() < min_num_hits_in_bundle) 
			{
				bb.clear();
				continue;
			}

			solve_bundle(bundle(bb), count, fout1, fout2);
			bb.clear();
			count++;
			if(max_num_bundles > 0 && count > max_num_bundles) break;
		}
		bb.add_hit(h, b);
    }

    bam_destroy1(b);
    bam_hdr_destroy(h);
    sam_close(fn);
	fout1.close();
	fout2.close();

	return 0;
}

int scallop::solve_bundle(const bundle &bd, int index, ofstream &fout1, ofstream &fout2)
{
	bd.print(index);

	splice_graph gr;
	bd.build_splice_graph(gr);

	assembler asmb(gr);

	asmb.solve();

	bd.output_gtf(fout1, asmb.paths0, index);
	bd.output_gtf(fout2, asmb.paths1, index);

	return 0;
}

int scallop::solve_gtf(const string &file)
{
	return 0;
}

int scallop::solve_example(const string &file)
{
	splice_graph gr;
	build_splice_graph(file, gr);
	assembler asmb(gr);
	asmb.draw("splice_graph.tex");
	asmb.solve();
	return 0;
}


