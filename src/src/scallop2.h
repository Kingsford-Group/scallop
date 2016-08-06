#ifndef __SCALLOP2_H__
#define __SCALLOP2_H__

#include "path.h"
#include "equation.h"
#include "splice_graph.h"
#include "nested_graph.h"
#include "hyper_edge.h"
#include "gateway.h"

typedef map< edge_descriptor, vector<int> > MEV;
typedef pair< edge_descriptor, vector<int> > PEV;
typedef pair< vector<int>, vector<int> > PVV;
typedef pair<PEE, int> PPEEI;
typedef map<PEE, int> MPEEI;
typedef pair<int, int> PI;
typedef map<int, int> MI;

// for noisy splice graph
class scallop2
{
public:
	scallop2();
	scallop2(const string &name, splice_graph &gr);
	scallop2(const string &name, splice_graph &gr, const vector<hyper_edge> &vhe);
	virtual ~scallop2();

public:
	string name;						// name for this gene
	splice_graph gr;					// splice graph

	vector<gateway> gateways;			// constructed gateway

	MEI e2i;							// edge map, from edge to index
	VE i2e;								// edge map, from index to edge
	MEV mev;							// super edges
	int round;							// round in iteration

	vector<path> paths;					// predicted transcripts

public:
	int clear();
	int save(scallop2 &sc);
	int load(scallop2 &sc);

public:
	int assemble();
	int assemble0();
	int assemble1();
	int assemble2();
	int greedy();

private:
	// trivial, or hard
	int classify();

	// init
	int init_super_edges();
	int init_gateways(const vector<hyper_edge> &vhe);

	// use hyper edges
	bool decompose_with_gateways();
	bool decompose_with_gateway(int v);

	// decompose trivial edges and vertices
	bool join_trivial_edges();
	bool join_trivial_edge(edge_descriptor &e);
	bool join_trivial_vertices();
	bool join_trivial_vertex(int i);
	bool decompose_trivial_vertices();
	bool decompose_trivial_vertex(int v);
	bool decompose_trivial_vertex(int v, vector<int> &ve);

	// compute shortest distances to source and target 
	int compute_shortest_source_distances();
	int compute_shortest_target_distances();
	int assign_reliability();

	// build equivalent classes
	bool infer_equivalent_classes();
	bool infer_equivalent_class(const vector<int> &v);
	bool infer_trivial_edges();
	bool infer_trivial_in_edge(int v);
	bool infer_trivial_out_edge(int v);
	bool infer_edges();
	bool infer_in_edges(int v);
	bool infer_out_edges(int v);
	bool infer_vertices();
	bool infer_vertex_with_in_edges(int v);
	bool infer_vertex_with_out_edges(int v);

	// rescale weights (coverage)
	int rescale_weights();
	int rescale_5end_weights(int i);
	int rescale_3end_weights(int i);

	// remove false boundary edges
	bool remove_false_boundary_edges();
	bool verify_false_boundary_edge(edge_descriptor e);

	// smooth weights
	bool smooth_splice_graph();
	bool smooth_splice_graph(const equation &eqn);
	bool smooth_splice_graph(const vector<equation> &eqns);
	bool smooth_vertex(int v);
	bool smooth_vertex(int v, const vector<equation> &eqns);
	bool smooth_trivial_equation(const equation &eqn);

	// iteratively decompose
	int iterate(bool greedy);
	bool decompose_with_equations(int level);

	// identify and handle equations 
	int identify_equations0(vector<equation> &eqns);
	int identify_equations1(vector<equation> &eqns);
	int identify_equations2(vector<equation> &eqns);
	int identify_equation(const vector<int> &subs, vector<equation> &eqns);
	bool verify_equation_nontrivial(equation &eqn);
	bool verify_equation_mergable(equation &eqn);

	// resolve equation
	int resolve_equation(equation &eqn);
	int resolve_equation(equation &eqn, vector<int> &ve);
	int resolve_equation(vector<int> &s, vector<int> &t, int &ma, int &md, vector<int> &ve);
	double compute_equation_error_ratio(equation &eqn);

	// split, and merge
	int split_edge(int exi, double w);
	int merge_adjacent_equal_edges(int x, int y);
	int merge_adjacent_edges(int x, int y);
	int split_merge_path(const vector<int> &p, double ww, vector<int> &vv);
	int split_merge_path(const VE &p, double ww, vector<int> &vv);

	// check, and make two edges adjacent
	bool check_adjacent_mergable(int ex, int ey, vector<PI> &p);
	bool check_adjacent_mergable(int ex, int ey, nested_graph &nt);
	int check_distant_mergable(int x, int y, double w, VE &p);
	int check_distant_mergable(int x, int y, double w);
	int build_adjacent_edges(const vector<PI> &p);

	// decompose the graph with greedy algorithm
	int greedy_decompose(int num);

	// collect existing s-t path e
	int collect_path(int e);
	int collect_existing_st_paths();

	// test, print and draw
	int draw_splice_graph(const string &file);
	int print();
};

#endif
