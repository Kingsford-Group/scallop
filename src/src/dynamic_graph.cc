#include "dynamic_graph.h"
#include "draw.h"
#include <sstream>
#include <fstream>

using namespace std;

using namespace dynamic_graph;

// basic operations
int dynamic_graph::add_vertex(splice_graph &gr) { return gr.add_vertex(); }
PEB dynamic_graph::add_edge(int s, int t, splice_graph &gr) { return PEB(gr.add_edge(s, t), true); }
PEB dynamic_graph::edge(int s, int t, const splice_graph &gr) { return gr.edge(s, t); }
int dynamic_graph::remove_edge(edge_descriptor e, splice_graph &gr) { return gr.remove_edge(e); }
int dynamic_graph::remove_edge(int s, int t, splice_graph &gr) { return gr.remove_edge(s, t); }
int dynamic_graph::clear_vertex(int v, splice_graph &gr) { return gr.clear_vertex(v); }
size_t dynamic_graph::num_vertices(const splice_graph &gr) { return gr.num_vertices(); }
size_t dynamic_graph::num_edges(const splice_graph &gr) { return gr.num_edges(); }
int dynamic_graph::source(edge_descriptor e, const splice_graph &gr) { return e->source(); }
int dynamic_graph::target(edge_descriptor e, const splice_graph &gr) { return e->target(); }
PEE dynamic_graph::in_edges(int v, const splice_graph &gr) { return gr.in_edges(v); }
PEE dynamic_graph::out_edges(int v, const splice_graph &gr) { return gr.out_edges(v); }
PEE dynamic_graph::edges(const splice_graph &gr) { return gr.edges(); }
int dynamic_graph::clear(splice_graph &gr) { return gr.clear(); }
int dynamic_graph::degree(int v, const splice_graph &gr) { return gr.degree(v); }
int dynamic_graph::in_degree(int v, const splice_graph &gr) { return gr.in_degree(v); }
int dynamic_graph::out_degree(int v, const splice_graph &gr) { return gr.out_degree(v); }
set<int> dynamic_graph::adjacent_vertices(int v, const splice_graph &gr) { return gr.adjacent_vertices(v); }
double dynamic_graph::get_vertex_weight(int v, const splice_graph &gr) { return gr.get_vertex_weight(v); }
double dynamic_graph::get_vertex_stddev(int v, const splice_graph &gr) { return gr.get_vertex_stddev(v); }
double dynamic_graph::get_edge_weight(edge_base *e, const splice_graph &gr) { return gr.get_edge_weight(e); }
double dynamic_graph::get_edge_stddev(edge_base *e, const splice_graph &gr) { return gr.get_edge_stddev(e); }
int dynamic_graph::set_vertex_weight(int v, double w, splice_graph &gr) { return gr.set_vertex_weight(v, w); }
int dynamic_graph::set_vertex_stddev(int v, double w, splice_graph &gr) { return gr.set_vertex_stddev(v, w); }
int dynamic_graph::set_edge_weight(edge_base *e, double w, splice_graph &gr) { return gr.set_edge_weight(e, w); }
int dynamic_graph::set_edge_stddev(edge_base *e, double w, splice_graph &gr) { return gr.set_edge_stddev(e, w); }

int dynamic_graph::build_splice_graph(splice_graph &gr, const string &file) { return gr.build(file); }
int dynamic_graph::write_splice_graph(const splice_graph &gr, const string &file) { return gr.write(file); }
int dynamic_graph::draw_splice_graph(const splice_graph &gr, const string &file, double len) { return gr.draw(file, len); }
int dynamic_graph::simulate_splice_graph(splice_graph &gr, int n, int m) { return gr.simulate(n, m); }
int dynamic_graph::compute_num_paths(const splice_graph &gr) { return gr.compute_num_paths(); }

int dynamic_graph::get_edge_weights(const splice_graph &gr, MED &med) { med = gr.get_edge_weights(); return 0; }
int dynamic_graph::set_edge_weights(splice_graph &gr, const MED &med) { return gr.set_edge_weights(med); }
int dynamic_graph::get_vertex_weights(const splice_graph &gr, vector<double> &v) { v = gr.get_vertex_weights(); return 0;}
int dynamic_graph::set_vertex_weights(splice_graph &gr, const vector<double> &v) { return gr.set_vertex_weights(v); }
int dynamic_graph::get_edge_indices(const splice_graph &gr, VE &i2e, MEI &e2i) { return gr.get_edge_indices(i2e, e2i); }
bool dynamic_graph::check_nested_splice_graph(const splice_graph &gr) { return gr.check_nested(); }
bool dynamic_graph::check_fully_connected(const splice_graph &gr) { return gr.check_fully_connected(); }
int dynamic_graph::bfs_splice_graph(const splice_graph &gr, int s, vector<int> &v) { return gr.bfs(s, v); }
