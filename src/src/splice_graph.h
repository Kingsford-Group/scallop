#ifndef __SPLICE_GRAPH_H__
#define __SPLICE_GRAPH_H__

// boost::graph
#include "boost/graph/graph_traits.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/properties.hpp"

using namespace boost;

//typedef property<vertex_index_t, int> vertex_properties;
//typedef property<edge_index_t, int> edge_properties;
//typedef adjacency_list<vecS, vecS, directedS> splice_graph;

//struct vertex_weight_t { typedef vertex_property_tag kind; };

namespace boost {
	enum vertex_weight_t { vertex_weight };
	enum vertex_stddev_t { vertex_stddev };
	enum edge_stddev_t { edge_stddev };
	//enum edge_weight_t { edge_weight };
	BOOST_INSTALL_PROPERTY(vertex, weight);
	BOOST_INSTALL_PROPERTY(vertex, stddev);
	BOOST_INSTALL_PROPERTY(edge, stddev);
	//BOOST_INSTALL_PROPERTY(edge, weight);
}

typedef property< vertex_stddev_t, double, property<vertex_weight_t, double> > vertex_properties;
typedef property< edge_stddev_t, double, property<edge_weight_t, double> > edge_properties;
typedef adjacency_list<vecS, vecS, bidirectionalS, vertex_properties, edge_properties> splice_graph;

typedef graph_traits<splice_graph>::vertex_iterator vertex_iterator;
typedef graph_traits<splice_graph>::vertex_descriptor vertex_descriptor;
typedef graph_traits<splice_graph>::in_edge_iterator in_edge_iterator;
typedef graph_traits<splice_graph>::out_edge_iterator out_edge_iterator;
typedef graph_traits<splice_graph>::edge_iterator edge_iterator;
typedef graph_traits<splice_graph>::edge_descriptor edge_descriptor;

typedef property_map<splice_graph, edge_weight_t>::type edge_weight_map;
//typedef property_map<splice_graph, vertex_index_t>::const_type const_vertex_index_map;

//static vertex_descriptor VNULL = graph_traits<splice_graph>::null_vertex();

#include <map>
using namespace std;

typedef map<edge_descriptor, bool> MEB;
typedef map<edge_descriptor, double> MED;
typedef pair<edge_descriptor, double> PED;
typedef pair<edge_descriptor, bool> PEB;
typedef map<edge_descriptor, int> MEI;
typedef pair<edge_descriptor, int> PEI;

int build_splice_graph(const string &file, splice_graph &gr);

#endif
