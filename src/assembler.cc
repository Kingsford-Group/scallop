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
#include "filter.h"
#include "cluster.h"

assembler::assembler()
{
	load_ccsread_info();
    sfn = sam_open(input_file.c_str(), "r");
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
		if((p.flag & 0x100) >= 1 && use_second_alignment == false) continue;	// secondary alignment
		if(p.n_cigar > MAX_NUM_CIGAR) continue;									// ignore hits with more than MAX_NUM_CIGAR cigar types
		if(p.qual < min_mapping_quality) continue;								// ignore hits with small quality
		if(p.n_cigar < 1) continue;												// should never happen

		hit ht(b1t);
		ht.set_tags(b1t);
		ht.set_strand();

		if(ccs.find(ht.qname) == ccs.end()) continue;

		ccsread_info cci = ccs[ht.qname];
		ht.set_ccsread_info(cci);

		//char rev = ((p.flag & 0x10) >= 1) ? '1' : '0';
		//printf("strand: rev = %c, ts = %c, xs = %c, ccstrand = %c\n", rev, ht.ts, ht.xs, cci.strand);

		// TODO, check rules for generating splice positions
		ht.build_splice_positions();

		//ht.print();
		//if(ht.nh >= 2 && p.qual < min_mapping_quality) continue;
		//if(ht.nm > max_edit_distance) continue;
		//if(uniquely_mapped_only == true && ht.nh != 1) continue; // TODO

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

		//printf("read strand = %c, xs = %c, ts = %c\n", ht.strand, ht.xs, ht.ts);

		// add hit
		if(ht.strand == '+') bb1.add_hit(ht);
		if(ht.strand == '-') bb2.add_hit(ht);
	}

	pool.push_back(bb1);
	pool.push_back(bb2);
	process(0);

	assign_RPKM();

	/*
	filter ft(trsts);
	ft.merge_single_exon_transcripts();
	trsts = ft.trs;
	*/

	write();
	
	return 0;
}

int assembler::process(int n)
{
	if(pool.size() < n) return 0;
	for(int i = 0; i < pool.size(); i++)
	{
		assemble(pool[i]);
		//exit(0);	// TODO
	}
	pool.clear();
	return 0;
}

int assembler::assemble(bundle_base &bb)
{
	if(bb.hits.size() < min_num_hits_in_bundle) return 0;
	if(bb.tid < 0) return 0;

	char buf[1024];
	strcpy(buf, hdr->target_name[bb.tid]);

	bundle bd(bb);

	bd.chrm = string(buf);
	bd.gr.gid = "gene." + tostring(index++);
	bd.build();

	bd.print(index);

	scallop sc(bd.gr, bd.hs);
	sc.assemble();

	if(verbose >= 2)
	{
		printf("transcripts:\n");
		for(int i = 0; i < sc.trsts.size(); i++) sc.trsts[i].write(cout);
	}

	cluster clst(sc.trsts);
	clst.solve();

	filter ft(clst.cct);
	ft.filter_length_coverage();
	//ft.remove_nested_transcripts();
	if(ft.trs.size() >= 1) trsts.insert(trsts.end(), ft.trs.begin(), ft.trs.end());

	if(verbose >= 2)
	{
		printf("transcripts after filtering (total transcripts = %lu):\n", trsts.size());
		for(int i = 0; i < ft.trs.size(); i++) ft.trs[i].write(cout);
	}

	printf("\n");
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

int assembler::load_ccsread_info()
{
	if(ccsread_info_file == "") return 0;
	ifstream fin(ccsread_info_file.c_str());
	if(fin.fail()) return 0;

	string s;
	while(getline(fin, s))
	{
		if(s.size() == 0) continue;
		ccsread_info cci(s);
		if(ccs.find(cci.qname) != ccs.end()) continue;
		ccs.insert(pair<string, ccsread_info>(cci.qname, cci));
	}

	/*
	map<string, ccsread_info>::iterator it;
	for(it = ccs.begin(); it != ccs.end(); it++)
	{
		it->second.print();
	}
	*/
	return 0;
}
