#ifndef __BUNDLE_H__
#define __BUNDLE_H__

#include "interval_map.h"
#include "bundle_base.h"
#include "junction.h"
#include "region.h"
#include "partial_exon.h"
#include "splice_graph.h"
#include "hyper_set.h"
#include "path.h"

using namespace std;

class bundle : public bundle_base
{
public:
	bundle(const bundle_base &bb);
	virtual ~bundle();

public:
	vector<junction> junctions;		// splice junctions
	vector<PPI> lsoft;				// left soft positions 
	vector<PPI> rsoft;				// right soft positions
	vector<PPI> lhard;				// left hard positions 
	vector<PPI> rhard;				// right hard positions
	splice_graph jr;				// junction graph
	vector<region> regions;			// regions
	vector<partial_exon> pexons;	// partial exons
	split_interval_map pmap;		// partial exon map
	splice_graph gr;				// splice graph
	hyper_set hs;					// hyper edges

public:
	virtual int build();
	int output_transcript(ofstream &fout, const path &p, const string &gid, const string &tid) const;	
	int output_transcripts(ofstream &fout, const vector<path> &p, const string &gid) const;	
	int print(int index);

private:
	// check whether hits are sorted
	int check_left_ascending();
	int check_right_ascending();

	// compute strand
	int compute_strand();

	// junction graph, for paired-end reads
	int build_junctions();
	int build_clips();
	int build_junction_graph();
	int draw_junction_graph(const string &file);
	int search_junction_graph(int32_t p);
	int traverse_junction_graph(int s, int t, VE &ve);
	int traverse_junction_graph1(int s, int t, VE &ve);
	int traverse_junction_graph1(int s, int t);
	int test_junction_graph();

	// build partial exons
	int align_hits();
	int build_regions();
	int build_partial_exons();
	int build_partial_exon_map();
	int locate_left_partial_exon(int32_t x);
	int locate_right_partial_exon(int32_t x);

	// super junctions and super partial_exons;
	int build_hyper_edges1();			// single end
	int build_hyper_edges2();			// paired end
	int build_splice_graph();			// weights from junctions
	int assign_edge_info_weights();		// weights from hyper edges
	int split_partial_exons();

	// store the corresponding pexons in each junction
	int link_partial_exons();
};


#endif
