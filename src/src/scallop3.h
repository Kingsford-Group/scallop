#ifndef __SCALLOP3_H__
#define __SCALLOP3_H__

#include "disjoint_sets.h"
#include "assembler.h"
#include "algebra.h"
#include "nested_graph.h"

typedef map< edge_descriptor, vector<int> > MEV;
typedef pair< edge_descriptor, vector<int> > PEV;
typedef pair<int, int> PI;

// algorithm: identify subsetsum signal
class scallop3 : public assembler
{
public:
	scallop3(const string &name, splice_graph &gr);
	virtual ~scallop3();

public:
	MEI e2i;				// edge map, from edge to index
	VE i2e;					// edge map, from index to edge
	MEV mev;				// super edges
	VVD ns;					// null space from nodes
	disjoint_sets_t ds;		// edges with the same weight are grouped together
	nested_graph nt;		// nested graph
	VE i2n;					// edge map, from index to edges in nt

public:
	int assemble();

private:
	// simplify the splice graph and init all data structures
	int init_super_edges();
	int reconstruct_splice_graph();
	bool decompose_trivial_vertex(int x);
	int init_disjoint_sets();
	int build_null_space();

	// get informations from ds since some edges are deleted
	vector<int> compute_representatives();
	vector< vector<int> > compute_disjoint_sets();

	// iteratively identify equations and update
	int iterate();
	int identify_equation(int &ei, vector<int> &sub);
	bool verify_equation(int ei, const vector<int> &sub);
	int split_edge(int ei, const vector<int> &sub);
	int build_nested_graph();
	bool identify_linkable_edges(int &ex, int &ey);
	bool connect_adjacent_edges(int x, int y);
	int compute_closest_subset(int xi, int w, const vector<PI> & xxp);
	int compute_closest_equal_edges(int &ex, int &ey);

	// test, print and draw
	int print();
	int draw_splice_graph(const string &file);
	int draw_nested_graph(const string &file);
};

#endif
