#include <cstdio>
#include <cassert>

#include "config.h"
#include "scallop.h"
#include "sgraph.h"

scallop::scallop()
{
}

scallop::~scallop()
{
}

int scallop::process(const string &file)
{
	string s = file.substr(file.size() - 3, 3);
	if(s == "bam" || s == "sam")
	{
		solve(file.c_str());
	}
	else
	{
		sgraph sg;
		sg.load(file);
		sg.draw("sgraph.tex");
		sg.solve();
	}
	return 0;
}

int scallop::solve_bundle(const bundle &bd, int index, ofstream &fout1, ofstream &fout2)
{
	bd.print(index);

	sgraph sg;
	sg.build(bd);
	//sg.draw("sgraph1.tex");

	sg.solve();

	//sg.draw("sgraph2.tex");
	bd.output_gtf(fout1, sg.paths0, index);
	bd.output_gtf(fout2, sg.paths1, index);

	return 0;
}

int scallop::solve(const char *bam_file)
{
    samFile *fn = sam_open(bam_file, "r");
    bam_hdr_t *h= sam_hdr_read(fn);
    bam1_t *b = bam_init1();

	ofstream fout1("greedy.gtf");
	ofstream fout2("iterat.gtf");

	int count = 0;
	bbase bb;
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
