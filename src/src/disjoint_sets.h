#ifndef __DISJOINT_SETS_H__
#define __DISJOINT_SETS_H__

#include <boost/pending/disjoint_sets.hpp>

using namespace boost;

typedef disjoint_sets_with_storage<identity_property_map, identity_property_map, find_with_path_halving> disjoint_sets_t;
//typedef disjoint_sets_with_storage<identity_property_map, identity_property_map, find_with_full_path_compression> ds_type;

int test_disjoint_sets();

#endif 
