#include <cstdio>
#include <cassert>
#include <sstream>

#include "config.h"
#include "assembler.h"
#include "bundle.h"
#include "scallop.h"
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
	else if(s == "gtf") assemble_gtf(input_file);
	else assemble_example(input_file);
	return 0;
}

int assembler::assemble_bam(const string &file)
{
    samFile *fn = sam_open(file.c_str(), "r");
    bam_hdr_t *h= sam_hdr_read(fn);
    bam1_t *b = bam_init1();

	ofstream fout;
	if(output_gtf_file != "") fout.open(output_gtf_file);

	int index = 0;
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

			bundle bd(bb);
			bd.print(index);
			
			splice_graph gr;
			bd.build_splice_graph(gr);

			string name = "bundle." + tostring(index);
			scallop sc(name, gr);
			sc.assemble();

			if(output_gtf_file != "")
			{
				bd.output_gtf(fout, sc.paths, algo, index);
			}

			index++;
			bb.clear();
			if(max_num_bundles > 0 && index > max_num_bundles) break;
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
		//if(name != "ZNF415") continue;

		splice_graph gr;
		gg.build_splice_graph(gr);

		string s;	
		int p0 = gr.compute_num_paths();
		int p1 = gr.num_edges() - gr.num_vertices() + 2;
		int p2 = gg.transcripts.size();
		assert(p0 >= p1);

		if(p0 == p1) assert(p2 >= p1);
		if(p0 == p1) s = "TRIVIAL";
		else if(p2 >= p1) s = "EASY";
		else s = "HARD";

		if(algo == "")
		{
			printf("gene %s, %lu transcipts, %d paths, %d required, %s\n", name.c_str(), gg.transcripts.size(), p0, p1, s.c_str());

			if(output_gtf_file == "") continue;
			if(s == "HARD") gg.output_gtf(fout);
		}

		if(s != "HARD") continue;

		scallop sc(name, gr);
		sc.assemble();

		if(output_gtf_file == "") continue;
		gg.output_gtf(fout, sc.paths, algo);
	}

	if(output_gtf_file != "") fout.close();

	return 0;
}

int assembler::assemble_example(const string &file)
{
	splice_graph gr;
	gr.build(file);

	scallop sc("example", gr);
	sc.assemble();

	return 0;
}
