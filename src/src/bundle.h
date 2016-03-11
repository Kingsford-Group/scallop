#ifndef __BUNDLE_H__
#define __BUNDLE_H__

#include "splice_graph.h"
#include "interval_map.h"
#include "bundle_base.h"
#include "junction.h"
#include "boundary.h"
#include "region.h"
#include "path.h"

using namespace std;
using namespace boost_splice_graph;

class bundle : public bundle_base
{
public:
	bundle(const bundle_base &bb);
	virtual ~bundle();

private:
	split_interval_map imap;				// interval map
	vector<junction> junctions;		// splice junctions
	vector<boundary> boundaries;	// all types of boundaries
	vector<region> regions;			// regions

public:
	int print(int index) const;
	int output_gtf(ofstream &fout, const vector<path> &paths, const string &prefix, int index) const;	
	int build_splice_graph(splice_graph &gr) const;

private:
	// check whether hits are sorted
	int check_left_ascending();
	int check_right_ascending();

	// build interval map
	int build_split_interval_map();

	// remove these intervals starting at a LEFT_BOUNDARY
	int remove_left_boundary_intervals();

	// binary search for a specific given starting point, return count
	int locate_hits(int32_t p, int &li);

	// infer boundaries
	int infer_junctions();
	int infer_left_boundaries();
	int infer_right_boundaries();
	int add_start_boundary();
	int add_end_boundary();

	// build regions
	int build_regions();
	// store the corresponding regions in each junction
	int link_regions();
	// run create_split on the boundaries of all regions
	int split_region_boundaries();
};

#endif
