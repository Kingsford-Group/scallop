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
	/*
	nested_graph nt;		// nested graph
	int gr_version;			// version of splice graph
	int nt_version;			// version of nested graph
	*/

	MEI e2i;				// edge map, from edge to index
	VE i2e;					// edge map, from index to edge
	MEV mev;				// super edges
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
	int smooth_with_equation(equation &eqn);

	// trivial, or hard
	int classify();

	// simplify the splice graph and init all data structures
	int init_super_edges();
	int reconstruct_splice_graph();
	bool init_trivial_vertex(int x);

	// remove empty edges, not use now
	int remove_empty_edges();

	// identify and handle equations 
	bool identify_equations1(vector<equation> &eqns);
	bool identify_equations2(vector<equation> &eqns);
	bool identify_equations3(vector<equation> &eqns);
	bool verify_equation_nontrivial(equation &eqn);
	bool verify_equation_mergable(equation &eqn);

	// resolve equation
	int resolve_equation(equation &eqn);
	int resolve_equation(vector<int> &s, vector<int> &t, int &ma, int &md);

	// use equation to decompose trivial vertex, new feature
	bool resolve_vertex_with_equations(vector<equation> &eqns);
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
