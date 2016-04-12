#include "bundle_extend.h"

bundle_extend::bundle_extend(const bundle_base &b, const genome &g)
	: gn(g), bundle(b)
{}

int bundle_extend::build()
{
	check_left_ascending();

	build_split_interval_map();

	infer_junctions();

	add_start_boundary();
	add_end_boundary();

	return 0;

	build_regions();
	link_regions();
	split_region_boundaries();

	return 0;
}

int bundle_extend::correct_boundaries()
{
	return 0;
}
