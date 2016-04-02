#include <cstdio>
#include <cassert>
#include <sstream>

#include "config.h"
#include "manager.h"
#include "bundle.h"
#include "stringtie.h"
#include "scallop.h"
#include "gtf_gene.h"
#include "nested_graph.h"

manager::manager()
{
}

manager::~manager()
{
}

int manager::process()
{
	if(input_file == "") return 0;
	string s = input_file.substr(input_file.size() - 3, 3);
	if(s == "bam" || s == "sam") assemble_bam(input_file);
	else if(s == "gtf") assemble_gtf(input_file);
	else assemble_example(input_file);
	return 0;
}

int manager::assemble_bam(const string &file)
{
    samFile *fn = sam_open(file.c_str(), "r");
    bam_hdr_t *h= sam_hdr_read(fn);
    bam1_t *b = bam_init1();

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

			/*
			bundle bd(bb);
			bd.print(index);
			
			splice_graph gr;
			bd.build_splice_graph(gr);

			stringtie st(gr);
			st.assemble();
			st.print("stringtie");

			scallop sc("", gr);
			sc.assemble();
			sc.print();

			bd.output_gtf(stringtie_fout, st.paths, "stringtie", index);
			bd.output_gtf(scallop_fout, sc.paths, "scallop", index);
			*/

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

int manager::assemble_gtf(const string &file)
{
	ifstream fin(file.c_str());
	if(fin.fail())
	{
		printf("open file %s error\n", file.c_str());
		return 0;
	}

	char line[102400];
	
	vector<gtf_gene> genes;
	map<string, int> m;
	while(fin.getline(line, 102400, '\n'))
	{
		gtf_exon ge(line);
		if(ge.feature != "exon") continue;
		if(m.find(ge.gene_id) == m.end())
		{
			gtf_gene gg;
			gg.add_exon(ge);
			genes.push_back(gg);
			m.insert(pair<string, int>(ge.gene_id, genes.size() - 1));
		}
		else
		{
			genes[m[ge.gene_id]].add_exon(ge);
		}
	}

	ofstream fout;
	if(output_gtf_file != "") fout.open(output_gtf_file);

	for(int i = 0; i < genes.size(); i++)
	{
		gtf_gene &gg = genes[i];
		if(gg.exons.size() <= 0) continue;

		string name = gg.exons[0].gene_id;

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
			printf("gene %s, %lu transcipts, total %lu exons, %d paths, %d required, %s\n", 
					name.c_str(), gg.transcripts.size(), gg.exons.size(), p0, p1, s.c_str());

			if(output_gtf_file == "") continue;
			if(s == "HARD") gg.output_gtf(fout);
			continue;
		}

		if(s != "HARD") continue;

		assembler * a;

		if(algo == "scallop") a = new scallop(name, gr);
		else if(algo == "stringtie") a = new stringtie(name, gr);
		else continue;

		a->assemble();

		if(output_gtf_file == "") continue;
		gg.output_gtf(fout, a->paths, algo);
	}

	if(output_gtf_file != "") fout.close();

	return 0;
}

int manager::assemble_example(const string &file)
{
	splice_graph gr;
	gr.build(file);

	scallop sc("example", gr);
	sc.assemble();

	return 0;
}
