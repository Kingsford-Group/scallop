#ifndef __SPLICE_GRAPH_H__
#define __SPLICE_GRAPH_H__

// boost::graph
#include "boost/graph/graph_traits.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/properties.hpp"
#include "boost/graph/adjacency_iterator.hpp"

using namespace boost;

#include <map>

using namespace std;

//typedef property<vertex_index_t, int> vertex_properties;
//typedef property<edge_index_t, int> edge_properties;
//typedef adjacency_list<vecS, vecS, directedS> splice_graph;

//struct vertex_weight_t { typedef vertex_property_tag kind; };

namespace boost 
{
	enum vertex_weight_t { vertex_weight };
	enum vertex_stddev_t { vertex_stddev };
	enum edge_stddev_t { edge_stddev };
	//enum edge_weight_t { edge_weight };
	BOOST_INSTALL_PROPERTY(vertex, weight);
	BOOST_INSTALL_PROPERTY(vertex, stddev);
	BOOST_INSTALL_PROPERTY(edge, stddev);
	//BOOST_INSTALL_PROPERTY(edge, weight);
}

namespace boost_graph
{
	typedef property< vertex_stddev_t, double, property<vertex_weight_t, double, property<vertex_color_t, default_color_type> > > vertex_properties;
	typedef property< edge_stddev_t, double, property<edge_weight_t, double> > edge_properties;
	typedef adjacency_list<multisetS, vecS, bidirectionalS, vertex_properties, edge_properties> splice_graph;

	typedef graph_traits<splice_graph>::vertex_iterator vertex_iterator;
	typedef graph_traits<splice_graph>::vertex_descriptor vertex_descriptor;
	typedef graph_traits<splice_graph>::in_edge_iterator in_edge_iterator;
	typedef graph_traits<splice_graph>::out_edge_iterator out_edge_iterator;
	typedef graph_traits<splice_graph>::edge_iterator edge_iterator;
	typedef graph_traits<splice_graph>::edge_descriptor edge_descriptor;
	typedef graph_traits<splice_graph>::adjacency_iterator adj_iterator;

	typedef property_map<splice_graph, edge_weight_t>::type edge_weight_map;
	//typedef property_map<splice_graph, vertex_index_t>::const_type const_vertex_index_map;

	//static vertex_descriptor VNULL = graph_traits<splice_graph>::null_vertex();

	typedef map<edge_descriptor, bool> MEB;
	typedef map<edge_descriptor, double> MED;
	typedef pair<edge_descriptor, double> PED;
	typedef pair<edge_descriptor, bool> PEB;
	typedef map<edge_descriptor, int> MEI;
	typedef pair<edge_descriptor, int> PEI;
	typedef vector<edge_descriptor> VE;

	// read, write, draw and simulate splice graph
	int build_splice_graph(splice_graph &gr, const string &file);
	int write_splice_graph(const splice_graph &gr, const string &file);
	int draw_splice_graph(const splice_graph &gr, const string &file, double len = 1.2);
	int simulate_splice_graph(splice_graph &gr, int n, int m);

	// get and set properties of splice graph
	int get_edge_weights(const splice_graph &gr, MED &med);
	int set_edge_weights(splice_graph &gr, const MED & med);
	int get_vertex_weights(const splice_graph &gr, vector<double> &v);
	int set_vertex_weights(splice_graph &gr, const vector<double> &v);
	int get_edge_indices(const splice_graph &gr, VE &i2e, MEI &e2i);

	// analysis the structure of splice graph
	int compute_num_paths(const splice_graph &gr);
	bool check_nested_splice_graph(const splice_graph &gr);
	bool check_directed_path(const splice_graph &gr, int s, int t);
	bool check_fully_connected(const splice_graph &gr);
	bool check_fully_reachable_from_start(const splice_graph &gr);
	bool check_fully_reachable_to_end(const splice_graph &gr);
	int bfs_distance(const splice_graph &gr, int s, vector<int> &v);


	// tests
	int test_bfs_distance();
	int test_remove_edge();
}

#endif
