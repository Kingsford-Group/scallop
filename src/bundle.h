/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

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
#include "gene.h"
#include "transcript.h"

using namespace std;

class bundle : public bundle_base
{
public:
	bundle(const bundle_base &bb);
	virtual ~bundle();

public:
	vector<junction> junctions;		// splice junctions
	vector<region> regions;			// regions
	vector<partial_exon> pexons;	// partial exons
	split_interval_map pmap;		// partial exon map
	splice_graph gr;				// splice graph
	hyper_set hs;					// hyper edges

public:
	virtual int build();
	int output_transcripts(ofstream &fout, const vector<path> &p, const string &gid) const;	
	int output_transcripts(gene &gn, const vector<path> &p, const string &gid) const;	
	int output_transcripts(vector<transcript> &trsts, const vector<path> &p, const string &gid) const;	
	int output_transcript(ofstream &fout, const path &p, const string &gid, const string &tid) const;	
	int output_transcript(transcript &trst, const path &p, const string &gid, const string &tid) const;	
	int count_junctions() const;
	int print(int index);

private:
	// check and init
	int check_left_ascending();
	int check_right_ascending();
	int compute_strand();

	// splice graph
	int build_junctions();
	int correct_junctions();
	int build_regions();
	int build_partial_exons();
	int link_partial_exons();
	int build_splice_graph();
	int build_partial_exon_map();
	int locate_left_partial_exon(int32_t x);
	int locate_right_partial_exon(int32_t x);

	// revise splice graph
	VE compute_maximal_edges();
	int revise_splice_graph();
	int refine_splice_graph();
	bool keep_surviving_edges();
	bool extend_boundaries();
	bool remove_small_junctions();
	bool remove_small_exons();
	bool remove_inner_boundaries();
	bool remove_intron_contamination();

	// super edges
	int build_hyper_edges2();			// paired end
	bool bridge_read(int x, int y, vector<int> &s);

};

#endif
