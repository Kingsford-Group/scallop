/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cassert>
#include <sstream>

#include "config.h"
#include "gtf.h"
#include "genome.h"
#include "assembler.h"
#include "scallop.h"
#include "sgraph_compare.h"
#include "super_graph.h"
#include "filter.h"

assembler::assembler(const config &c)
  : cfg(c)
{
    sfn = sam_open(cfg.input_file.c_str(), "r");
    hdr = sam_hdr_read(sfn);
    b1t = bam_init1();
	index = 0;
	terminate = false;
	qlen = 0;
	qcnt = 0;
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

		if((p.flag & 0x4) >= 1) continue;										// read is not mapped
		if((p.flag & 0x100) >= 1 && cfg.use_second_alignment == false) continue;	// secondary alignment
		if(p.n_cigar > MAX_NUM_CIGAR) continue;									// ignore hits with more than 7 cigar types
		if(p.qual < cfg.min_mapping_quality) continue;								// ignore hits with small quality
		if(p.n_cigar < 1) continue;												// should never happen

		hit ht(b1t, cfg);
		ht.set_tags(b1t);
		ht.set_strand();
		ht.build_splice_positions();

		//ht.print();

		//if(ht.nh >= 2 && p.qual < min_mapping_quality) continue;
		//if(ht.nm > max_edit_distance) continue;
		//if(ht.verify_junctions() == false) continue;

		qlen += ht.qlen;
		qcnt += 1;

		// truncate
		if(ht.tid != bb1.tid || ht.pos > bb1.rpos + cfg.min_bundle_gap)
		{
			pool.push_back(bb1);
			bb1.clear();
		}
		if(ht.tid != bb2.tid || ht.pos > bb2.rpos + cfg.min_bundle_gap)
		{
			pool.push_back(bb2);
			bb2.clear();
		}

		// process
		process(cfg.batch_bundle_size);

		// add hit
		if(cfg.uniquely_mapped_only == true && ht.nh != 1) continue;
		if(cfg.library_type != UNSTRANDED && ht.strand == '+' && ht.xs == '-') continue;
		if(cfg.library_type != UNSTRANDED && ht.strand == '-' && ht.xs == '+') continue;
		if(cfg.library_type != UNSTRANDED && ht.strand == '.' && ht.xs != '.') ht.strand = ht.xs;
		if(cfg.library_type != UNSTRANDED && ht.strand == '+') bb1.add_hit(ht);
		if(cfg.library_type != UNSTRANDED && ht.strand == '-') bb2.add_hit(ht);
		if(cfg.library_type == UNSTRANDED && ht.xs == '.') bb1.add_hit(ht);
		if(cfg.library_type == UNSTRANDED && ht.xs == '.') bb2.add_hit(ht);
		if(cfg.library_type == UNSTRANDED && ht.xs == '+') bb1.add_hit(ht);
		if(cfg.library_type == UNSTRANDED && ht.xs == '-') bb2.add_hit(ht);
	}

	pool.push_back(bb1);
	pool.push_back(bb2);
	process(0);

	assign_RPKM();

	filter ft(trsts, cfg);
	ft.merge_single_exon_transcripts();
	trsts = ft.trs;

	write();

	return 0;
}

int assembler::process(int n)
{
	if(pool.size() < n) return 0;

	for(int i = 0; i < pool.size(); i++)
	{
		bundle_base &bb = pool[i];
		if(bb.hits.size() < cfg.min_num_hits_in_bundle) continue;
		if(bb.tid < 0) continue;

		char buf[1024];
		strcpy(buf, hdr->target_name[bb.tid]);

		bundle bd(bb, cfg);

		bd.chrm = string(buf);
		bd.build();
		if(cfg.verbose >= 1) bd.print(index);

    parameter_advising<bundle, transcript_set, config_list>(bd, cfg.gene_config_list);
    assemble(bd.gr, bd.hs);
		index++;
	}
	pool.clear();
	return 0;
}

int assembler::assemble(const splice_graph &gr0, const hyper_set &hs0)
{
	super_graph sg(gr0, hs0, cfg);
	sg.build();

	vector<transcript> gv;
	for(int k = 0; k < sg.subs.size(); k++)
	{
		string gid = "gene." + tostring(index) + "." + tostring(k);
		if(cfg.fixed_gene_name != "" && gid != cfg.fixed_gene_name) continue;

		if(cfg.verbose >= 2 && (k == 0 || cfg.fixed_gene_name != "")) sg.print();

		splice_graph &gr = sg.subs[k];
		hyper_set &hs = sg.hss[k];

		gr.gid = gid;
		scallop sc(gr, hs, cfg);
		sc.assemble();

		if(cfg.verbose >= 2)
		{
			printf("transcripts:\n");
			for(int i = 0; i < sc.trsts.size(); i++) sc.trsts[i].write(cout);
		}

		filter ft(sc.trsts, cfg);
		ft.join_single_exon_transcripts();
		ft.filter_length_coverage();
		if(ft.trs.size() >= 1) gv.insert(gv.end(), ft.trs.begin(), ft.trs.end());

		if(cfg.verbose >= 2)
		{
			printf("transcripts after filtering:\n");
			for(int i = 0; i < ft.trs.size(); i++) ft.trs[i].write(cout);
		}

		if(cfg.fixed_gene_name != "" && gid == cfg.fixed_gene_name) terminate = true;
		if(terminate == true) return 0;
	}

	filter ft(gv, cfg);
	ft.remove_nested_transcripts();
	if(ft.trs.size() >= 1) trsts.insert(trsts.end(), ft.trs.begin(), ft.trs.end());

	return 0;
}

int assembler::assign_RPKM()
{
	double factor = 1e9 / qlen;
	for(int i = 0; i < trsts.size(); i++)
	{
		trsts[i].assign_RPKM(factor);
	}
	return 0;
}

int assembler::write()
{
	ofstream fout(cfg.output_file.c_str());
	if(fout.fail()) return 0;
	for(int i = 0; i < trsts.size(); i++)
	{
		transcript &t = trsts[i];
		t.write(fout);
	}
	fout.close();
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
