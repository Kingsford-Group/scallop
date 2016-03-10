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
	vector<bool> i2b;		// indicate whether this edge has been deleted
	MEV mev;				// super edges
	VVD ns;					// null space from nodes
	disjoint_sets_t ds;		// edges with the same weight are grouped together

public:
	int assemble();
	int print();

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
	int split_edge(int ei, const vector<int> &sub);
	bool connect_adjacent_edges(int x, int y);
	int compute_closest_subset(int xi, int w, const vector<PI> & xxp);
	int compute_closest_equal_edges(int &ex, int &ey);
};

#endif
