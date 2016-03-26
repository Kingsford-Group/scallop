#ifndef __SCALLOP_H__
#define __SCALLOP_H__

#include "disjoint_sets.h"
#include "assembler.h"
#include "nested_graph.h"

typedef map< edge_descriptor, vector<int> > MEV;
typedef pair< edge_descriptor, vector<int> > PEV;
typedef pair<int, int> PI;

// algorithm: identify subsetsum signal
class scallop : public assembler
{
public:
	scallop(const string &name, splice_graph &gr);
	virtual ~scallop();

public:
	MEI e2i;				// edge map, from edge to index
	VE i2e;					// edge map, from index to edge
	MEV mev;				// super edges
	disjoint_sets_t ds;		// edges with the same weight are grouped together
	nested_graph nt;		// nested graph
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
	set<int> compute_singletons();

	// main flow of the algorithm
	bool iterate();

	// identify equations 
	bool identify_equation(int &ei, vector<int> &sub);
	bool identify_edge_equation(int ei, vector<int> &sub);

	// split exi w.r.t eyi
	int split_edge(int exi, double w);
	int split_edge(int exi, int eyi);
	vector<int> split_edge(int ei, const vector<int> &sub);

	// check, and make two equal edges adjacent 
	bool identify_linkable_edges(int &ex, int &ey, vector<PI> &p);
	bool check_linkable(int ex, int ey, vector<PI> &p);
	int build_adjacent_equal_edges(const vector<PI> &p);
	int connect_adjacent_equal_edges(int x, int y);

	// compute, and merge two closest equal edges with minimum distance
	bool compute_closest_equal_edges(int &ex, int &ey);
	int connect_equal_edges(int x, int y);

	// decompose trivial vertices 
	bool decompose_trivial_vertices();

	// greedily decompose the splice graph into paths
	double compute_maximum_weight_path(vector<int> &p);

	// split the path out
	int connect_path(const vector<int> &p, double ww);

	// test, print and draw
	int draw_splice_graph(const string &file);
	int print();
};

#endif
