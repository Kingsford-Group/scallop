/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cassert>
#include <sstream>

#include "gtf.h"
#include "genome1.h"
#include "config.h"
#include "assembler.h"
#include "scallop.h"
#include "sgraph_compare.h"
#include "super_graph.h"
#include "filter.h"
#include "estimator.h"

assembler::assembler()
{
    sfn = sam_open(input_file.c_str(), "r");
    hdr = sam_hdr_read(sfn);
    b1t = bam_init1();
	terminate = false;
	index = -1;
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
		if((p.flag & 0x100) >= 1 && use_second_alignment == false) continue;	// secondary alignment
		if(p.n_cigar > MAX_NUM_CIGAR) continue;									// ignore hits with more than 7 cigar types
		if(p.qual < min_mapping_quality) continue;								// ignore hits with small quality
		if(p.n_cigar < 1) continue;												// should never happen

		hit ht(b1t);
		//if(ht.nm > max_edit_distance) continue;
		//if(ht.verify_junctions() == false) continue;

		qlen += ht.qlen;
		qcnt += 1;

		// truncate
		if(ht.tid != bb1.tid || ht.pos > bb1.rpos + min_bundle_gap)
		{
			pool.push_back(bb1);
			bb1.clear();
		}
		if(ht.tid != bb2.tid || ht.pos > bb2.rpos + min_bundle_gap)
		{
			pool.push_back(bb2);
			bb2.clear();
		}

		// process
		process(batch_bundle_size);

		// add hit
		if(uniquely_mapped_only == true && ht.nh != 1) continue;
		if(library_type != UNSTRANDED && ht.strand == '.') continue;
		if(library_type != UNSTRANDED && ht.strand == '+' && ht.xs == '-') continue;
		if(library_type != UNSTRANDED && ht.strand == '-' && ht.xs == '+') continue;
		if(library_type != UNSTRANDED && ht.strand == '+') bb1.add_hit(ht);
		if(library_type != UNSTRANDED && ht.strand == '-') bb2.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '.') bb1.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '.') bb2.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '+') bb1.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '-') bb2.add_hit(ht);
	}

	pool.push_back(bb1);
	pool.push_back(bb2);
	process(0);

	//merge_multi_exon_transcripts();
	assign_RPKM();
	write();
	
	return 0;
}

int assembler::process(int n)
{
	if(pool.size() < n) return 0;
	for(int i = 0; i < pool.size(); i++)
	{
		assemble(pool[i]);
	}
	pool.clear();
	return 0;
}

int assembler::assemble(const bundle_base &bb)
{
	//printf("hits.size = %lu, min-num = %d\n", bb.hits.size(), min_num_hits_in_bundle);

	if(bb.hits.size() < min_num_hits_in_bundle) return 0;
	if(bb.tid < 0) return 0;

	char buf[1024];
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

	/*
	if(ref_file != "") compare(bd.gr, ref_file, "compare.tex");
	if(ref_file1 != "" && bd.strand == '+') compare(bd.gr, ref_file1, "compare1.tex");
	if(ref_file2 != "" && bd.strand == '-') compare(bd.gr, ref_file2, "compare2.tex");
	*/

	for(int k = 0; k < sg.subs.size(); k++)
	{
		//if(k != maxk) continue;

		string gid = "bundle." + tostring(index) + "." + tostring(k);
		if(fixed_gene_name != "" && gid != fixed_gene_name) continue;

		if(verbose >= 1 && (k == 0 || fixed_gene_name != "")) bd.print(index);
		if(verbose >= 2 && (k == 0 || fixed_gene_name != "")) sg.print();

		if(algo == "empty") continue;

		splice_graph &gr = sg.subs[k];
		hyper_set &hs = sg.hss[k];

		/*
		if(ref_file != "") compare(gr, ref_file, "compare.tex");
		if(ref_file1 != "" && bd.strand == '+') compare(gr, ref_file1, "compare1.tex");
		if(ref_file2 != "" && bd.strand == '-') compare(gr, ref_file2, "compare2.tex");
		*/

		//if(gr.num_vertices() <= 3 && sg.subs.size() >= 2) continue;
		//if(gr.num_vertices() <= 3 && bd.junctions.size() >= 1) continue;

		scallop sc(gid, gr, hs);
		sc.assemble();

		/*
		estimator est(gr, sc.paths);
		est.estimate();
		*/

		vector<path> pp;
		for(int i = 0; i < sc.paths.size(); i++)
		{
			path p;
			p.v = sg.get_root_vertices(k, sc.paths[i].v);
			p.abd = sc.paths[i].abd;
			p.reads = sc.paths[i].reads;
			pp.push_back(p);
		}

		vector<transcript> vv;
		bd.output_transcripts(vv, pp, gid);

		filter ft(vv);
		ft.join();
		ft.select();
		if(ft.trs.size() >= 1) trsts.insert(trsts.end(), ft.trs.begin(), ft.trs.end());

		if(fixed_gene_name != "" && gid == fixed_gene_name) terminate = true;
		if(terminate == true) return 0;
	}

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

int assembler::merge_multi_exon_transcripts()
{
	vector<transcript> v1;
	vector<transcript> v2;
	for(int i = 0; i < trsts.size(); i++)
	{
		if(trsts[i].exons.size() <= 1) v1.push_back(trsts[i]);
		else v2.push_back(trsts[i]);
	}

	genome1 gm;
	gm.build(v2);
	gm.merge();

	for(int i = 0; i < gm.transcripts.size(); i++)
	{
		v1.push_back(gm.transcripts[i]);
	}

	trsts = v1;
	return 0;
}

int assembler::write()
{
	ofstream fout(output_file.c_str());
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
