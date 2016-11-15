#ifndef __SCALLOP3_H__
#define __SCALLOP3_H__

#include "splice_graph.h"
#include "hyper_set.h"
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
typedef pair<int, double> PID;
typedef map<int, double> MID;
typedef pair<PI, double> PPID;

// for noisy splice graph
class scallop
{
public:
	scallop();
	scallop(const string &name, const splice_graph &gr, const hyper_set &hs);
	virtual ~scallop();

public:
	int assemble();

public:
	string name;						// name for this gene
	splice_graph gr;					// splice graph
	MEI e2i;							// edge map, from edge to index
	VE i2e;								// edge map, from index to edge
	MEV mev;							// super edges
	MED med;							// reads for this super edge
	vector<int> v2v;					// vertex map
	hyper_set hs;						// hyper edges
	int round;							// round in iteration
	vector<path> paths;					// predicted transcripts
	vector<router> routers;

private:
	// init
	int classify();
	int init_vertex_map();
	int init_super_edges();
	int init_inner_weights();
	int add_pseudo_hyper_edges();
	int assert_weights();
	int refine_splice_graph();

	// resolve iteratively
	bool resolve_small_edges();
	bool resolve_hyper_tree(int status);
	bool resolve_hyper_vertex(int status);
	bool resolve_trivial_vertex();
	bool resolve_hyper_edge1();
	bool resolve_hyper_edge0();

	// smooth vertex
	int balance_vertex(int x);
	double compute_balance_ratio(int x);
	double balance_vertex(undirected_graph &ug, const vector<int> &u2e, vector<PPID> &vpi);
	int complete_graph(undirected_graph &ug, const vector<int> &u2e, int root);

	int get_weights(int v, MID &m);
	int set_weights(MID &m);

	// decomposing subroutines
	int decompose_tree(const vector<PPID> &vpi);
	int decompose_trivial_vertex(int v);
	int split_vertex(int x, const vector<int> &xe, const vector<int> &ye);
	int split_edge(int exi, double w);
	int merge_adjacent_edges(int x, int y);
	int merge_adjacent_edges(int x, int y, double ww);
	int merge_adjacent_equal_edges(int x, int y);
	int remove_edge(int e);
	int split_merge_path(const VE &p, double w);
	int split_merge_path(const vector<int> &p, double w);
	int collect_path(int e);
	int collect_existing_st_paths();
	int greedy_decompose(int num);
	double compute_smallest_edge(int x, int &e);
	double compute_smallest_removable_edge(int &se);
	double compute_smallest_splitable_vertex(int &root, int status);

	// print and draw
	int print();
	int stats();
	int draw_splice_graph(const string &file);
	vector<int> topological_sort();
};

#endif
