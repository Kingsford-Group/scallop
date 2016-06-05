#ifndef __SCALLOP2_H__
#define __SCALLOP2_H__

#include "splice_graph.h"
#include "nested_graph.h"
#include "path.h"
#include "equation.h"

typedef map< edge_descriptor, vector<int> > MEV;
typedef pair< edge_descriptor, vector<int> > PEV;
typedef pair< vector<int>, vector<int> > PVV;
typedef pair<int, int> PI;

// for perfectly estimated splice graph
class scallop
{
public:
	scallop();
	scallop(const string &name, splice_graph &gr);
	virtual ~scallop();

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
	int clear();
	int save(scallop &sc);
	int load(scallop &sc);

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
	bool decompose_with_equations();
	int decompose_single_equation(equation &eqn);

	// trivial, or hard
	int classify();

	// simplify the splice graph and init all data structures
	int init_super_edges();
	int reconstruct_splice_graph();
	bool init_trivial_vertex(int x);

	// remove empty edges
	int remove_empty_edges();

	// identify and handle equations 
	int identify_equation1(vector<equation> &eqns);
	int identify_equation(const vector<int> &subs, vector<int> &subt);
	bool verify_equation_nontrivial(const vector<int> &subs, const vector<int> &subt);

	// split and connect pairs of identical edges
	PI connect_pairs(const vector<int> &vx, const vector<int> &vy);

	// split exi w.r.t eyi
	int split_edge(int exi, double w);
	int split_edge(int exi, int eyi);
	vector<int> split_edge(int ei, const vector<int> &sub);

	// check, and make two adjacent equal edges adjacent 
	bool check_adjacent_mergable(int ex, int ey, vector<PI> &p);
	bool check_adjacent_mergable(int ex, int ey);
	int build_adjacent_equal_edges(const vector<PI> &p);
	int connect_adjacent_equal_edges(int x, int y);

	// compute, and merge two distant equal edges
	int check_distant_mergable(int x, int y, VE &p);
	int check_distant_mergable(int x, int y);
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
