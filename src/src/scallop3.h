#ifndef __SCALLOP3_H__
#define __SCALLOP3_H__

#include "splice_graph.h"
#include "equation.h"
#include "router.h"
#include "path.h"

typedef map< edge_descriptor, vector<int> > MEV;
typedef pair< edge_descriptor, vector<int> > PEV;
typedef pair< vector<int>, vector<int> > PVV;
typedef pair<PEE, int> PPEEI;
typedef map<PEE, int> MPEEI;
typedef pair<int, int> PI;
typedef map<int, int> MI;

// for noisy splice graph
class scallop3
{
public:
	scallop3();
	scallop3(const string &name, const splice_graph &gr);
	virtual ~scallop3();

public:
	string name;						// name for this gene
	splice_graph gr;					// splice graph

	vector<router> routers;			// constructed router

	MEI e2i;							// edge map, from edge to index
	VE i2e;								// edge map, from index to edge
	MEV mev;							// super edges
	int round;							// round in iteration

	vector<path> paths;					// predicted transcripts

public:
	int assemble();

private:
	// trivial, or hard
	int classify();

	// decompose
	int iterate();
	bool decompose();
	bool decompose_trivial_vertex(int v);
	bool decompose_trivial_vertex(int v, vector<int> &ve);
	int split_edge(int exi, double w);
	int merge_adjacent_edges(int x, int y);
	int merge_adjacent_equal_edges(int x, int y);

	// smooth vertex
	bool smooth_vertex(int v, const vector<equation> &eqns);

	// init
	int init_super_edges();
	int init_routers(const vector<hyper_edge> &vhe);

	// print and draw
	int print();
	int draw_splice_graph(const string &file);
};

#endif
