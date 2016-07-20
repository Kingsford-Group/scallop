#ifndef __SCALLOP2_H__
#define __SCALLOP2_H__

#include "path.h"
#include "equation.h"
#include "splice_graph.h"
#include "nested_graph.h"

typedef map< edge_descriptor, vector<int> > MEV;
typedef pair< edge_descriptor, vector<int> > PEV;
typedef pair< vector<int>, vector<int> > PVV;
typedef pair<int, int> PI;
typedef map<int, int> MI;

// for noisy splice graph
class scallop2
{
public:
	scallop2();
	scallop2(const string &name, splice_graph &gr);
	virtual ~scallop2();

public:
	string name;			// name for this gene
	splice_graph gr;		// splice graph

	MEI e2i;				// edge map, from edge to index
	VE i2e;					// edge map, from index to edge
	MEV mev;				// super edges
	int round;				// round in iteration

	vector<path> paths;		// predicted transcripts

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

	// decompose trivial edges and vertices
	bool decompose_trivial_edges();
	bool decompose_trivial_vertices();
	bool decompose_trivial_vertex(int i);
	bool decompose_trivial_edge(edge_descriptor &e);

	// smooth weights
	int smooth();
	int smooth(equation &eqn);

	// iteratively decompose
	int iterate();
	bool decompose_with_equations(int level);

	// identify and handle equations 
	int identify_equations1(vector<equation> &eqns);
	int identify_equations2(vector<equation> &eqns);
	int identify_equation(const vector<int> &subs, vector<equation> &eqns);
	bool verify_equation_nontrivial(equation &eqn);
	bool verify_equation_mergable(equation &eqn);

	// resolve equation
	int resolve_equation(equation &eqn);
	int resolve_equation(vector<int> &s, vector<int> &t, int &ma, int &md);

	// use equation to decompose trivial vertex, new feature
	bool resolve_vertex_with_equations(vector<equation> &eqns);
	bool resolve_vertex_with_equation(equation &eqn);
	bool resolve_vertex_with_equation1(equation &eqn);
	bool resolve_vertex_with_equation2(equation &eqn);

	// split, and merge
	int split_edge(int exi, double w);
	int merge_adjacent_equal_edges(int x, int y);
	int split_merge_path(const vector<int> &p, double ww, vector<int> &vv);
	int split_merge_path(const VE &p, double ww, vector<int> &vv);
	int merge_adjacent_edges(int x, int y);

	// check, and make two edges adjacent
	bool check_adjacent_mergable(int ex, int ey, vector<PI> &p);
	bool check_adjacent_mergable(int ex, int ey, nested_graph &nt);
	int check_distant_mergable(int x, int y, double w, VE &p);
	int check_distant_mergable(int x, int y, double w);
	int build_adjacent_edges(const vector<PI> &p);

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
