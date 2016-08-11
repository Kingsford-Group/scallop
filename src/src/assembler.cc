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
#include "sgraph_compare.h"

assembler::assembler()
{
}

assembler::~assembler()
{
}

int assembler::process()
{
	if(input_file == "")
	{
		if(output_file == "") return 0;
		splice_graph sg;
		sg.simulate(simulation_num_vertices, simulation_num_edges, simulation_max_edge_weight);
		sg.write(output_file);
		return 0;
	}

	string s = input_file.substr(input_file.size() - 3, 3);
	if(s == "bam" || s == "sam") assemble_bam(input_file);
	if(s == "gtf") assemble_gtf(input_file);
	if(s == "sgr") assemble_sgr(input_file);
	return 0;
}

int assembler::assemble_bam(const string &file)
{
    samFile *fn = sam_open(file.c_str(), "r");
    bam_hdr_t *h= sam_hdr_read(fn);
    bam1_t *b = bam_init1();

	ofstream fout;
	if(output_file != "")
	{
		fout.open(output_file);
		if(fout.fail()) return 0;
	}

	int index = -1;
	bundle_base bb1;		// for + reads
	bundle_base bb2;		// for - reads

	bool terminate = false;
    while(sam_read1(fn, h, b) >= 0)
	{
		bam1_core_t &p = b->core;

		if((p.flag & 0x4) >= 1) continue;			// read is not mapped
		if((p.flag & 0x100) >= 1) continue;			// secondary alignment
		if(p.n_cigar < 1) continue;					// should never happen
		if(p.n_cigar > MAX_NUM_CIGAR) continue;		// ignore hits with more than 7 cigar types
		//if(p.qual <= min_quality_score) continue;	// ignore hits with quality-score < 5
		
		if(bb1.get_num_hits() > 0 && (bb1.get_rpos() + min_bundle_gap < p.pos || p.tid != bb1.get_tid()))
		{
			int f = process_bundle(bb1, h, index, fout);
			if(f == 1) break;
		}

		if(bb2.get_num_hits() > 0 && (bb2.get_rpos() + min_bundle_gap < p.pos || p.tid != bb2.get_tid()))
		{
			int f = process_bundle(bb2, h, index, fout);
			if(f == 1) break;
		}

		hit ht(b);

		if(ht.strand == '+') bb1.add_hit(ht);
		if(ht.strand == '-') bb2.add_hit(ht);
    }
	process_bundle(bb1, h, index, fout);
	process_bundle(bb2, h, index, fout);

    bam_destroy1(b);
    bam_hdr_destroy(h);
    sam_close(fn);
	if(output_file != "") fout.close();

	return 0;
}

int assembler::process_bundle(bundle_base &bb, bam_hdr_t *h, int &index, ofstream &fout)
{
	if(bb.get_num_hits() < min_num_hits_in_bundle) 
	{
		bb.clear();
		return 0;
	}

	char buf[1024];
	int tid = bb.get_tid();
	assert(tid >= 0);
	strcpy(buf, h->target_name[tid]);
	bb.set_chrm(buf);

	index++;

	string name = "bundle." + tostring(index);

	if(fixed_gene_name != "" && name != fixed_gene_name)
	{
		bb.clear();
		return 0;
	}

	bundle bd(bb);
	bd.build();

	bb.clear();

	bd.print(index);

	if(bd.size() >= 100) return 0;

	splice_graph gr;
	vector<hyper_edge> vhe;
	bd.build_splice_graph(gr, vhe);

	scallop2 sc(name, gr, vhe);
	sc.assemble();

	if(ref_file != "")
	{
		genome g(ref_file);
		gtf gg(g.genes[0]);
		splice_graph gt;
		gg.build_splice_graph(gt);

		sgraph_compare sgc(gt, gr);
		sgc.compare();
	}

	if(output_file != "") bd.output_gtf(fout, sc.paths, algo, index);

	if(fixed_gene_name != "" && name == fixed_gene_name) return 1;
	else return 0;
}

int assembler::assemble_gtf(const string &file)
{
	genome g(file);

	ofstream fout;
	if(output_file != "") fout.open(output_file);

	for(int i = 0; i < g.genes.size(); i++)
	{
		gtf gg(g.genes[i]);

		if(gg.transcripts.size() < min_gtf_transcripts_num) continue;

		string name = gg.get_gene_id();

		// DEBUG
		if(fixed_gene_name != "" && name != fixed_gene_name) continue;

		splice_graph gr;
		gg.build_splice_graph(gr);

		scallop1 sc(name, gr);
		sc.assemble();

		if(output_file == "") continue;
		gg.output_gtf(fout, sc.paths, algo);
	}

	if(output_file != "") fout.close();

	return 0;
}

int assembler::assemble_sgr(const string &file)
{
	splice_graph sg;
	sg.build(file);

	scallop1 sc("shao", sg);
	sc.assemble();

	return 0;
}
