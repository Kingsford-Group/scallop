#ifndef __SCALLOP3_H__
#define __SCALLOP3_H__

#include "disjoint_sets.h"
#include "assembler.h"
#include "algebra.h"

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
	set<PI> sis;			// set of intersecting edges
	disjoint_sets_t ds;		// edges with the same weight are grouped together
	int round;				// round in iteration

public:
	int assemble();

private:
	// simplify the splice graph and init all data structures
	int init_super_edges();
	int reconstruct_splice_graph();
	bool init_trivial_vertex(int x);
	int init_disjoint_sets();

	// get informations from ds since some edges are deleted
	vector<int> compute_representatives();
	vector< vector<int> > compute_disjoint_sets();

	// iteratively identify equations and update
	bool iterate();
	bool identify_equation(int &ei, vector<int> &sub);
	bool verify_equation(int ei, const vector<int> &sub);
	int split_edge(int exi, int eyi);
	vector<int> split_edge(int ei, const vector<int> &sub);
	int compute_intersecting_edges();
	bool identify_linkable_edges(int &ex, int &ey, vector<int> &p);
	bool check_linkable(int ex, int ey, vector<int> &p);
	int build_adjacent_edges(int ex, int ey, const vector<int> &p);
	int connect_adjacent_edges(int x, int y);
	bool decompose_trivial_vertices();
	bool compute_shortest_equal_edges(int &ex, int &ey);
	int connect_equal_edges(int x, int y);

	// test, print and draw
	int print();
	int draw_splice_graph(const string &file);
};

#endif
