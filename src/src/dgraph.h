#ifndef __DGRAPH_H__
#define __DGRAPH_H__

// boost::graph
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>

using namespace boost;

//typedef property<vertex_index_t, int> vertex_properties;
//typedef property<edge_index_t, int> edge_properties;
//typedef adjacency_list<vecS, vecS, directedS> dgraph;
typedef adjacency_list<vecS, vecS, bidirectionalS> dgraph;

typedef graph_traits<dgraph>::vertex_iterator vertex_iterator;
typedef graph_traits<dgraph>::vertex_descriptor vertex_descriptor;
typedef graph_traits<dgraph>::edge_iterator edge_iterator;
typedef graph_traits<dgraph>::out_edge_iterator out_edge_iterator;
typedef graph_traits<dgraph>::edge_descriptor edge_descriptor;

//typedef property_map<dgraph, vertex_index_t>::const_type const_vertex_index_map;
//typedef property_map<dgraph, edge_index_t>::type edge_index_map;

//static vertex_descriptor VNULL = graph_traits<dgraph>::null_vertex();

#include <map>
using namespace std;

typedef map<edge_descriptor, double> MED;
typedef map<edge_descriptor, int> MEI;
typedef pair<edge_descriptor, int> PEI;

#endif
