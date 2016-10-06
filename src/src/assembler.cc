#include <cstdio>
#include <cassert>
#include <sstream>

#include "genome.h"
#include "gtf.h"
#include "config.h"
#include "assembler.h"
#include "scallop3.h"
#include "sgraph_compare.h"
#include "super_graph.h"

assembler::assembler()
{
    sfn = sam_open(input_file.c_str(), "r");
    hdr = sam_hdr_read(sfn);
    b1t = bam_init1();
	terminate = false;
	index = -1;
	qlen = 0;
}

assembler::~assembler()
{
    bam_destroy1(b1t);
    bam_hdr_destroy(hdr);
    sam_close(sfn);
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
		qlen += ht.qlen;
		if(ht.strand == '.') continue;	// TODO

		truncate(ht);
		add_hit(ht);
    }

	for(int i = 0; i < vbb.size(); i++) process(vbb[i]);

	if(output_file == "") return 0;

	double factor = 1e9 * average_read_length / qlen;
	gm.assign_RPKM(factor);
	gm.write(output_file);

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
		if(ht.pos <= bb.rpos + min_bundle_gap && ht.tid == bb.tid) continue;

		process(bb);
		vbb.erase(vbb.begin() + i);
		truncate(ht);
	}
	return 0;
}

int assembler::process(const bundle_base &bb)
{
	printf("hits.size = %lu, min-num = %d\n", bb.hits.size(), min_num_hits_in_bundle);

	if(bb.hits.size() < min_num_hits_in_bundle) return 0;

	char buf[1024];
	assert(bb.tid >= 0);
	strcpy(buf, hdr->target_name[bb.tid]);

	bundle bd(bb);

	bd.chrm = string(buf);
	bd.build();

	index++;

	/*
	bd.print(index);
	if(ref_file != "") compare(bd.gr, ref_file, "compare.tex");
	if(ref_file1 != "" && bd.strand == '+') compare(bd.gr, ref_file1, "compare1.tex");
	if(ref_file2 != "" && bd.strand == '-') compare(bd.gr, ref_file2, "compare2.tex");
	*/

	super_graph sg(bd.gr, bd.hs);
	sg.build();

	int maxk = -1;
	int maxv = 0;
	for(int k = 0; k < sg.subs.size(); k++)
	{
		if(sg.subs[k].num_vertices() < maxv) continue;
		maxk = k;
		maxv = sg.subs[k].num_vertices();
	}

	if(ref_file != "") compare(bd.gr, ref_file, "compare.tex");
	if(ref_file1 != "" && bd.strand == '+') compare(bd.gr, ref_file1, "compare1.tex");
	if(ref_file2 != "" && bd.strand == '-') compare(bd.gr, ref_file2, "compare2.tex");

	for(int k = 0; k < sg.subs.size(); k++)
	{
		//if(k != maxk) continue;

		string gid = "bundle." + tostring(index) + "." + tostring(k);
		if(fixed_gene_name != "" && gid != fixed_gene_name) continue;

		if(k == 0 || fixed_gene_name != "") bd.print(index);
		if(k == 0 || fixed_gene_name != "") sg.print();

		splice_graph &gr = sg.subs[k];
		hyper_set &hs = sg.hss[k];

		/*
		if(ref_file != "") compare(gr, ref_file, "compare.tex");
		if(ref_file1 != "" && bd.strand == '+') compare(gr, ref_file1, "compare1.tex");
		if(ref_file2 != "" && bd.strand == '-') compare(gr, ref_file2, "compare2.tex");
		*/

		double ss = gr.compute_coverage() / average_read_length;
		if(ss < min_splice_graph_coverage) continue;

		scallop3 sc(gid, gr, hs);
		sc.assemble();

		vector<path> pp;
		for(int i = 0; i < sc.paths.size(); i++)
		{
			path p;
			p.v = sg.get_root_vertices(k, sc.paths[i].v);
			p.abd = sc.paths[i].abd;
			p.reads = sc.paths[i].reads;
			pp.push_back(p);
		}

		gene gn;
		bd.output_transcripts(gn, pp, gid);
		gm.add_gene(gn);

		if(fixed_gene_name != "" && gid == fixed_gene_name) terminate = true;
		if(terminate == true) return 0;
	}

	return 0;
}

int assembler::compare(splice_graph &gr, const string &file, const string &texfile)
{
	if(file == "") return 0;

	genome g(file);
	if(g.genes.size() <= 0) return 0;

	gtf gg(g.genes[0]);

	splice_graph gt;
	gg.build_splice_graph(gt);

	sgraph_compare sgc(gt, gr);
	sgc.compare(texfile);

	return 0;
}
