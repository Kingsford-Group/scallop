#include <cstdio>
#include <cassert>
#include <sstream>

#include "config.h"
#include "assembler.h"
#include "bundle.h"
#include "scallop1.h"
#include "scallop2.h"
#include "gtf.h"
#include "genome.h"
#include "nested_graph.h"

assembler::assembler()
{
}

assembler::~assembler()
{
}

int assembler::process()
{
	if(input_file == "") return 0;
	string s = input_file.substr(input_file.size() - 3, 3);
	if(s == "bam" || s == "sam") assemble_bam(input_file);
	if(s == "gtf") assemble_gtf(input_file);
	return 0;
}

int assembler::assemble_bam(const string &file)
{
    samFile *fn = sam_open(file.c_str(), "r");
    bam_hdr_t *h= sam_hdr_read(fn);
    bam1_t *b = bam_init1();

	ofstream fout;
	if(output_gtf_file != "") fout.open(output_gtf_file);

	int index = -1;
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

			index++;

			// DEBUG
			/*
			if(index <= 9000) 
			{
				bb.clear();
				continue;
			}
			*/

			bundle bd(bb);
			bd.build();

			bb.clear();
			if(max_num_bundles > 0 && index > max_num_bundles) break;

			bd.print(index);

			if(bd.size() >= 100) continue;

			splice_graph gr;
			bd.build_splice_graph(gr);

			string name = "bundle." + tostring(index);
			scallop2 sc(name, gr);
			sc.assemble();

			if(output_gtf_file != "") bd.output_gtf(fout, sc.paths, algo, index);

		}
		bb.add_hit(h, b);
    }

    bam_destroy1(b);
    bam_hdr_destroy(h);
    sam_close(fn);

	return 0;
}

int assembler::assemble_gtf(const string &file)
{
	genome g(file);

	ofstream fout;
	if(output_gtf_file != "") fout.open(output_gtf_file);

	for(int i = 0; i < g.genes.size(); i++)
	{
		gtf gg(g.genes[i]);
		if(gg.transcripts.size() <= 0) continue;

		string name = gg.get_gene_id();

		// DEBUG
		if(fixed_gene_name != "" && name != fixed_gene_name) continue;

		splice_graph gr;
		gg.build_splice_graph(gr);

		scallop1 sc(name, gr);
		sc.assemble();

		if(output_gtf_file == "") continue;
		gg.output_gtf(fout, sc.paths, algo);
	}

	if(output_gtf_file != "") fout.close();

	return 0;
}

