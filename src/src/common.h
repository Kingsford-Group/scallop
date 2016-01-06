#ifndef __COMMON_H__
#define __COMMON_H__

#include <map>
#include <stdint.h>
using namespace std;

// macros: using int64_t for two int32_t
#define pack(x, y) (int64_t)((((int64_t)(x)) << 32) | ((int64_t)(y)))
#define high32(x) (int32_t)((x) >> 32)
#define low32(x) (int32_t)(((x) << 32) >> 32)


// definitions
typedef map<int32_t, int> MPI;
typedef pair<int32_t, int> PPI;


// boost::interval map
#include <boost/icl/interval_map.hpp>
using namespace boost;

typedef pair<size_t, size_t> PT;
typedef icl::right_open_interval<int32_t> ROI;
typedef icl::interval_map<int32_t, int32_t, icl::partial_absorber, less, icl::inplace_plus, icl::inter_section, ROI> imap_t;

// return the overlap at position p
int compute_overlap(const imap_t &imap, int32_t p);
// return the sum of overlap in [p,q)
int cumulate_overlap(const imap_t &imap, int32_t p, int32_t q, int32_t t);


// boost::binomial distribution
#include <boost/math/distributions/binomial.hpp>
using namespace boost::math;

// return the score(transformed from probability)
// that >= x is observed
uint32_t compute_binomial_score(int n, double pr, int x);


// boost::graph
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>

//typedef property<vertex_index_t, int> vertex_properties;
//typedef property<edge_index_t, int> edge_properties;
typedef adjacency_list<vecS, vecS, directedS> dgraph;

typedef graph_traits<dgraph>::vertex_iterator vertex_iterator;
typedef graph_traits<dgraph>::vertex_descriptor vertex_descriptor;
typedef graph_traits<dgraph>::edge_iterator edge_iterator;
typedef graph_traits<dgraph>::out_edge_iterator out_edge_iterator;
typedef graph_traits<dgraph>::edge_descriptor edge_descriptor;

//typedef property_map<dgraph, vertex_index_t>::const_type const_vertex_index_map;
//typedef property_map<dgraph, edge_index_t>::type edge_index_map;

//static vertex_descriptor VNULL = graph_traits<dgraph>::null_vertex();

typedef map<edge_descriptor, double> MED;

#endif
