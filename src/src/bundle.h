#ifndef __BUNDLE_H__
#define __BUNDLE_H__

#include "imap.h"
#include "bbase.h"
#include "bridge.h"
#include "boundary.h"
#include "region.h"

using namespace std;

class bundle : public bbase
{
public:
	bundle(const bbase &bb);
	virtual ~bundle();

protected:
	imap_t imap;					// interval map
	vector<bridge> bridges;			// splice bridges
	vector<boundary> boundaries;	// all types of boundaries
	vector<region> regions;			// regions

private:
	// check whether hits are sorted
	int check_left_ascending();
	int check_right_ascending();

	// build interval map
	int build_interval_map();

	// remove these intervals starting at a LEFT_BOUNDARY
	int remove_left_boundary_intervals();

	// binary search for a specific given starting point, return count
	int locate_hits(int32_t p, int &li);

	// infer boundaries
	int infer_bridges();
	int infer_left_boundaries();
	int infer_right_boundaries();
	int add_start_boundary();
	int add_end_boundary();

	// build regions
	int build_regions();

	// store the corresponding regions in each bridge
	int link_regions();
};

#endif
