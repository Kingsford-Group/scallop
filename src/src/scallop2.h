#ifndef __SCALLOP2_H__
#define __SCALLOP2_H__

#include "splice_graph.h"
#include "nested_graph.h"
#include "path.h"

typedef map< edge_descriptor, vector<int> > MEV;
typedef pair< edge_descriptor, vector<int> > PEV;
typedef pair<int, int> PI;

// for perfectly estimated splice graph
class scallop2
{
public:
	scallop2(const string &name, splice_graph &gr);
	virtual ~scallop2();

public:
	string name;			// name for this gene
	splice_graph gr;		// splice graph

	MEI e2i;				// edge map, from edge to index
	VE i2e;					// edge map, from index to edge
	MEV mev;				// super edges
	nested_graph nt;		// nested graph
	int round;				// round in iteration

	vector<path> paths;		// predicted transcripts

public:
	int assemble();
	int assemble0();
	int assemble1();
	int assemble2();
	int greedy();

private:
	// iteratively decompose
	int iterate();
	bool decompose_trivial_vertices();
	bool decompose_with_equation();

	// trivial, or hard
	int classify();

	// simplify the splice graph and init all data structures
	int init_super_edges();
	int reconstruct_splice_graph();
	bool init_trivial_vertex(int x);

	// identify and handle equations 
	bool identify_equation1(vector<int> &subs, vector<int> &subt);
	int identify_equation(const vector<int> &subs, vector<int> &subt);
	bool verify_equation_nontrivial(const vector<int> &subs, const vector<int> &subt);

	// split and connect pairs of identical edges
	int connect_pairs(const vector<int> &vx, const vector<int> &vy);

	// split exi w.r.t eyi
	int split_edge(int exi, double w);
	int split_edge(int exi, int eyi);
	vector<int> split_edge(int ei, const vector<int> &sub);

	// check, and make two equal edges adjacent 
	bool check_linkable(int ex, int ey, vector<PI> &p);
	int build_adjacent_equal_edges(const vector<PI> &p);
	int connect_adjacent_equal_edges(int x, int y);

	// compute, and merge two closest equal edges with minimum distance
	bool connect_equal_edges(int x, int y);
	int connect_path(const vector<int> &p, double ww);
	int connect_path(const VE &p, double ww);

	// decompose the graph with greedy algorithm
	int greedy_decompose();

	// collect existing s-t path e
	int collect_path(int e);
	int collect_existing_st_paths();

	// test, print and draw
	int draw_splice_graph(const string &file);
	int print();
};

#endif
