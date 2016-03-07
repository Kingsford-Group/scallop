#ifndef __DISJOINT_SETS_H__
#define __DISJOINT_SETS_H__

#include <vector>
#include <boost/pending/disjoint_sets.hpp>

using namespace std;
using namespace boost;

typedef disjoint_sets_with_storage<identity_property_map, identity_property_map, find_with_path_halving> disjoint_sets_t;
//typedef disjoint_sets_with_storage<identity_property_map, identity_property_map, find_with_full_path_compression> ds_type;

// process functions
vector<int> get_representatives(disjoint_sets_t &ds, int n);
vector< vector<int> > get_disjoint_sets(disjoint_sets_t &ds, int n);

// test function
int test_disjoint_sets();

#endif 
