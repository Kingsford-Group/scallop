#ifndef __BUNDLE_H__
#define __BUNDLE_H__

#include "super_graph.h"
#include "interval_map.h"
#include "bundle_base.h"
#include "junction.h"
#include "partial_exon.h"
#include "hyper_edge.h"
#include "path.h"
#include "super_region.h"

using namespace std;

class bundle : public bundle_base
{
public:
	bundle(const bundle_base &bb);
	virtual ~bundle();

public:
	vector<junction> junctions;		// splice junctions
	splice_graph jr;				// junction graph
	vector<region> regions;			// regions
	vector<partial_exon> pexons;	// partial exons
	super_graph sg;					// super_graph

public:
	virtual int build();
	int output_transcript(ofstream &fout, const path &p, const string &gid, const string &tid) const;	
	int print(int index);

private:
	// check whether hits are sorted
	int check_left_ascending();
	int check_right_ascending();

	// compute strand
	int compute_strand();

	// junction graph, for paired-end reads
	int build_junctions();
	int build_junction_graph();
	int draw_junction_graph(const string &file);
	int search_junction_graph(int32_t p);
	int traverse_junction_graph(int s, int t, VE &ve);
	int traverse_junction_graph1(int s, int t, VE &ve);
	int traverse_junction_graph1(int s, int t);
	int test_junction_graph();

	// process hits
	int align_hits();
	bool verify_unique_mapping(const hit &h);
	int compute_read1_intervals(const hit &h, vector<int64_t> &vv);

	// build partial exons
	int build_regions();
	int build_partial_exons();

	// super junctions and super partial_exons;
	int build_hyper_edges(vector<hyper_edge> &vhe);
	int build_hyper_edges0(vector<hyper_edge> &vhe);
	int build_splice_graph(splice_graph &gr);
	int search_partial_exons(int32_t p);

	// store the corresponding pexons in each junction
	int link_partial_exons();
};


#endif
