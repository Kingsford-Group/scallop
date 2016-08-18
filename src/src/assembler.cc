#include <cstdio>
#include <cassert>
#include <sstream>

#include "genome.h"
#include "gtf.h"
#include "config.h"
#include "assembler.h"
#include "scallop2.h"
#include "sgraph_compare.h"

assembler::assembler()
{
    sfn = sam_open(input_file.c_str(), "r");
    hdr = sam_hdr_read(sfn);
    b1t = bam_init1();
	terminate = false;
	index = -1;
	if(output_file != "") fout.open(output_file);
}

assembler::~assembler()
{
    bam_destroy1(b1t);
    bam_hdr_destroy(hdr);
    sam_close(sfn);
	if(output_file != "") fout.close();
}

int assembler::assemble()
{
    while(sam_read1(sfn, hdr, b1t) >= 0)
	{
		if(terminate == true) return 0;

		bam1_core_t &p = b1t->core;

		if((p.flag & 0x4) >= 1) continue;			// read is not mapped
		if((p.flag & 0x100) >= 1) continue;		// secondary alignment
		if(p.n_cigar < 1) continue;				// should never happen
		if(p.n_cigar > MAX_NUM_CIGAR) continue;	// ignore hits with more than 7 cigar types
		if(p.qual < min_mapping_quality) continue;	// ignore hits with small quality

		hit ht(b1t);
		if(ht.strand == '.') continue;	// TODO

		truncate(ht);
		add_hit(ht);
    }

	for(int i = 0; i < vbb.size(); i++) process(vbb[i]);

	return 0;
}

int assembler::add_hit(const hit &ht)
{
	bool b = false;
	for(int i = 0; i < vbb.size(); i++)
	{
		bundle_base &bb = vbb[i];
		assert(bb.hits.size() >= 1);
		assert(bb.strand != '.');

		//if(bb.overlap(ht) == false) continue;
		if(ht.strand != bb.strand) continue;
		if(ht.tid != bb.tid) continue;
		if(ht.pos > bb.rpos + min_bundle_gap) continue;

		bb.add_hit(ht);
		b = true;
		break;
	}

	if(b == true) return 0;

	bundle_base bb;
	bb.add_hit(ht);

	vbb.push_back(bb);

	return 0;
}

int assembler::truncate(const hit &ht)
{
	for(int i = 0; i < vbb.size(); i++)
	{
		bundle_base &bb = vbb[i];
		//if(ht.pos < bb.rpos && ht.tid == bb.tid) continue;
		if(ht.pos <= bb.rpos + min_bundle_gap && ht.tid == bb.tid) continue;

		process(bb);
		vbb.erase(vbb.begin() + i);
		truncate(ht);
	}
	return 0;
}

int assembler::process(const bundle_base &bb)
{
	if(bb.hits.size() < min_num_hits_in_bundle) return 0;

	char buf[1024];
	assert(bb.tid >= 0);
	strcpy(buf, hdr->target_name[bb.tid]);

	bundle bd(bb);

	bd.chrm = string(buf);
	bd.build();

	if(bd.junctions.size() <= 0 && ignore_single_exon_transcripts) return 0;

	index++;
	bd.print(index);

	for(int k = 0; k < bd.sg.subs.size(); k++)
	{
		string name = "bundle." + tostring(index) + "." + tostring(k);

		if(fixed_gene_name != "" && name != fixed_gene_name) continue;

		splice_graph &gr = bd.sg.subs[k];

		if(ref_file != "") compare(gr);

		scallop2 sc(name, gr);
		sc.assemble();

		if(output_file != "")
		{
			for(int i = 0; i < sc.paths.size(); i++)
			{
				string tid = name + "." + tostring(i);
				path p;
				p.v = bd.sg.get_root_vertices(k, sc.paths[i].v);
				p.abd = sc.paths[i].abd;
				bd.output_transcript(fout, p, name, tid);
			}
		}

		if(fixed_gene_name != "" && name == fixed_gene_name) terminate = true;

		if(terminate == true) return 0;
	}

	return 0;
}

int assembler::compare(splice_graph &gr)
{
	if(ref_file == "") return 0;

	genome g(ref_file);
	if(g.genes.size() <= 0) return 0;

	gtf gg(g.genes[0]);

	splice_graph gt;
	gg.build_splice_graph(gt);

	sgraph_compare sgc(gt, gr);
	sgc.compare("compare.tex");
	//sgc.compare(string("compare.") + tostring(index) + string(".tex"));

	return 0;
}
