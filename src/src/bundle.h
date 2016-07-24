#ifndef __BUNDLE_H__
#define __BUNDLE_H__

#include "splice_graph.h"
#include "interval_map.h"
#include "bundle_base.h"
#include "junction.h"
#include "region.h"
#include "partial_exon.h"
#include "path.h"

using namespace std;

class bundle : public bundle_base
{
public:
	bundle(const bundle_base &bb);
	virtual ~bundle();

private:
	vector<junction> junctions;		// splice junctions
	splice_graph jr;				// junction graph
	map<int, int> p2v;				// position to vertex of jr
	split_interval_map imap;		// interval map
	vector<region> regions;			// regions
	vector<partial_exon> pexons;	// partial exons

public:
	virtual int build();
	int build_splice_graph(splice_graph &gr) const;
	int output_gtf(ofstream &fout, const vector<path> &paths, const string &prefix, int index) const;	
	int print(int index) const;
	int size() const;

protected:
	// check whether hits are sorted
	int check_left_ascending();
	int check_right_ascending();

	// infer boundaries
	int infer_junctions();

	// junction graph, for paired-end reads
	int build_junction_graph();
	int draw_junction_graph(const string &file);

	// build interval map
	int build_split_interval_map();

	// binary search for a specific given starting point, return count
	int locate_hits(int32_t p, int &li);

	// build regions
	int build_partial_exons();

	// store the corresponding regions in each junction
	int link_partial_exons();

	// run create_split on the boundaries of all regions
	// not necessary, because now using split_interval_map
	//int split_boundaries();
};

#endif
