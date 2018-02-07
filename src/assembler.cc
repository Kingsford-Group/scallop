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
  mapped_counts.push_back(0);
  mapped_counts.push_back(0);
  mapped_counts.push_back(0);
  mapped_counts.push_back(0);
  mapped_counts.push_back(0);
  mapped_counts.push_back(0);
  mapped_counts.push_back(0);
  mapped_counts.push_back(0);
  mapped_counts.push_back(0);
}

assembler::~assembler()
{
    bam_destroy1(b1t);
    bam_hdr_destroy(hdr);
    sam_close(sfn);
}

assembler* assembler::solve(const config &c){
  assembler* a = new assembler(c);
  a->assemble();
  return a;
}

bool assembler::operator>(const assembler &a){
  cout << mapped_counts[0] << "\t" << a.mapped_counts[0] << endl;
  return (mapped_counts[0] > a.mapped_counts[0]);
}

int assembler::assemble()
{
  uint32_t max_hit_index = 0;
  //vector<uint32_t> hit_index;
  //hit_index.reserve(1000000);
  unordered_map <string, bool> read_mapping; //(max_hit_index,-1);
  int bundle_index = -1;
  vector<hit> hits;
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

    //uint32_t hit_index = stoul(ht.qname.substr(ht.qname.find(".")+1));
    //if(hit_index > max_hit_index) max_hit_index = hit_index;
    //hit_index.push_back(stoul(ht.qname.substr(ht.qname.find(".")+1)));
    //if(hit_index[hit_index.size()-1] > max_hit_index) max_hit_index = hit_index[hit_index.size()-1];

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
    //read_mapping.resize(max_hit_index+1,-1);
		process(cfg.batch_bundle_size, bundle_index, read_mapping);

		// add hit
		bool add_hit = false;

    if(cfg.uniquely_mapped_only == true && ht.nh != 1) continue;
		if(cfg.library_type != UNSTRANDED && ht.strand == '+' && ht.xs == '-') continue;
		if(cfg.library_type != UNSTRANDED && ht.strand == '-' && ht.xs == '+') continue;
		if(cfg.library_type != UNSTRANDED && ht.strand == '.' && ht.xs != '.') ht.strand = ht.xs;
		if(cfg.library_type != UNSTRANDED && ht.strand == '+'){
      bb1.add_hit(ht);
      add_hit = true;
    }if(cfg.library_type != UNSTRANDED && ht.strand == '-'){
      bb2.add_hit(ht);
      add_hit = true;
    }if(cfg.library_type == UNSTRANDED && ht.xs == '.'){
      bb1.add_hit(ht);
      add_hit = true;
    }if(cfg.library_type == UNSTRANDED && ht.xs == '.'){
      bb2.add_hit(ht);
      add_hit = true;
    }if(cfg.library_type == UNSTRANDED && ht.xs == '+'){
      bb1.add_hit(ht);
      add_hit = true;
    }if(cfg.library_type == UNSTRANDED && ht.xs == '-'){
      bb2.add_hit(ht);
      add_hit = true;
    }

    if(add_hit){
      //mapped_counts[4] ++;
      //hits.push_back(ht);
      //if(bb1.tid > 0) hits[hits.size()-1].seqname = hdr->target_name[bb1.tid];
      //else if(bb2.tid > 0) hits[hits.size()-1].seqname = hdr->target_name[bb2.tid];
      //cerr << hits[hits.size()-1].qname << "\thit tid: " << hits[hits.size()-1].tid << "\thit seqname: " << hits[hits.size()-1].seqname << "\tbb tid: " << bb1.tid << "," << bb2.tid << endl;
    }
  }


	pool.push_back(bb1);
	pool.push_back(bb2);
	process(0, bundle_index, read_mapping);

	assign_RPKM();

	filter ft(trsts, cfg);
	ft.merge_single_exon_transcripts();
	trsts = ft.trs;
  for(int j=0; j<trsts.size(); j++){
    trsts[j].used = true;
    //trsts[j].gene_id += "-";
    //trsts[j].gene_id += to_string(j);
    //trsts[j].transcript_id += "-";
    //trsts[j].transcript_id += to_string(j);
  }
  //count_mapped(hits);
  cout << "\tmapped reads\total_edit_distance\tavg_edit_distance\ttotal_norm_edit_distance\tavg_norm_edit_distance\tpairs_mapped\ttotal_reads\ttotal_pairs\ttotal_variance\ttotal_variance_from_predicted\ttotal_bundles\n";
  cout << "Mapped counts:\t" << mapped_counts[0] << "\t" << mapped_counts[1] << "\t" << (1.0 * mapped_counts[1]/mapped_counts[0]) << "\t" << mapped_counts[2] << "\t" << (1.0 * mapped_counts[2]/mapped_counts[0]) << "\t" << mapped_counts[3] << "\t" << mapped_counts[4] << "\t" << mapped_counts[5] << "\t" << (mapped_counts[6]/1000.0) << "\t" << (mapped_counts[7]/1000.0) << "\t" << mapped_counts[8] <<  endl;
	return 0;
}

void assembler::count_mapped(vector<hit> hits){
	int mapped = 0;
  float total_distance = 0;
  float norm_distance = 0;
  int mapped_pairs = 0;
  int pairs = 0;
  int reads = 0;
  unordered_map<string,int> hit_index;
  unordered_map <string, int> pair_mapping; //(max_hit_index,-1);
  int bundle_index = -1;
  mapped_counts[4] = hits.size();
  for(int i=0; i<hits.size(); i++){
    cerr << "SEQNAME: " << hits[i].seqname << "\ttid: " << hits[i].tid << endl;
    string name = hits[i].qname;
    name += ((hits[i].flag & 0x40)?".1":".2");
    assert(!((hits[i].flag & 0x40) && (hits[i].flag & 0x80)));
    assert((hits[i].flag & 0x40) || (hits[i].flag & 0x80));

    bool found = false;
    bool has_pair = (hit_index.find(hits[i].qname) != hit_index.end());
    int ip = hit_index[hits[i].qname];

    for(int j=0; j<trsts.size(); j++){
      if(!trsts[j].used) continue;
			if(hits[i].maps_to_transcript(trsts[j])){
        cerr << name << "\t" << bundle_index << "\t" << trsts[j].seqname << "\t" << trsts[j].gene_id << "\t" << trsts[j].transcript_id << endl;
        if(!found){
          //cerr << name << "\t" << trsts[j].seqname << "\t" << trsts[j].gene_id << "\t" << trsts[j].transcript_id<< endl;
          mapped_counts[0]++;
          mapped_counts[1] += hits[i].nm;
          mapped_counts[2] += 1.0 * hits[i].nm / hits[i].qlen;
        }
        //pair_mapping[name] = 1;
        found = true;
				//mapped++;
				//break;

        if(has_pair){
          //if(pair_mapping.find(hits[i].qname) == pair_mapping.end()){
            if(hits[ip].maps_to_transcript(trsts[j])){
              mapped_counts[3]++;
              //pair_mapping[hits[i].qname] = 1;
              break;
            }
          //}
        }
      }
		}
    if(!found) cerr << name << "\t" << bundle_index << "\t-\t-\t-\n";
    if(has_pair) mapped_counts[5]++;
    /*string seen_name = "seen";
    seen_name += name;
    if(pair_mapping.find(seen_name) == pair_mapping.end()){
      mapped_counts[4]++;
      if(!found) cerr << name << "\t" << bundle_index << "\t-\t-\t-\n";
    }
    pair_mapping[seen_name] = 1;
    if(has_pair){
      seen_name = "seen";
      seen_name += hits[i].qname;
      if(pair_mapping.find(seen_name) == pair_mapping.end())
        mapped_counts[5]++;
      pair_mapping[seen_name] = 1;
    }*/
    hit_index[hits[i].qname] = i;
  }
}

int assembler::process(int n, int bundle_index, unordered_map <string, bool> &read_mapping)
{
  if(pool.size() < n) return 0;

  //cout << "Read Mapping Size: " << read_mapping.size() << endl;

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

    assemble(bd.gr, bd.hs, bb.trsts);
    bundle_index++;
    vector<uint64_t> counts = bb.hits_on_transcripts(bundle_index, read_mapping);
    bb.trsts.clear();

    //cout << "Done Clear" << endl;
    mapped_counts[0] += counts[0];
    mapped_counts[1] += counts[1];
    mapped_counts[2] += counts[2];
    mapped_counts[3] += counts[3];
    mapped_counts[4] += counts[4];
    mapped_counts[5] += counts[5];
    if(counts[6] != 9999999999){
      //cout << "counts[6]:" << counts[6] << endl;
      mapped_counts[6] += counts[6];
      mapped_counts[7] += counts[7];
      mapped_counts[8] += 1;
    }
    //cout << "Done assign" << endl;
    
		index++;
	}
	pool.clear();
	return 0;
}

int assembler::assemble(const splice_graph &gr0, const hyper_set &hs0, vector<transcript> &rtn)
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
	if(ft.trs.size() >= 1){
    trsts.insert(trsts.end(), ft.trs.begin(), ft.trs.end());
    rtn.insert(rtn.end(), ft.trs.begin(), ft.trs.end());
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

int assembler::write(const char* fname)
{
  cout << "Got here with output file: " << fname << endl;
	//ofstream fout(cfg.output_file.c_str());
	ofstream fout(fname);
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
